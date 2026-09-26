"""Host-to-host adaptive Lorenz solves through public Python solver APIs."""
import os
import time

os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")

import diffrax
import gradsolve
import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)


class Lorenz:
    name = "diffeqgpu_lorenz"
    dim = 3
    t0 = 0.0
    t1 = 1.0
    is_stiff = False

    @staticmethod
    def f_jax(t, y, p):
        return jnp.stack((10 * (y[1] - y[0]),
                          p[0] * y[0] - y[1] - y[0] * y[2],
                          y[0] * y[1] - (8 / 3) * y[2]))


def prepare(rhos, precision, rtol, library, device="gpu"):
    dtype = np.dtype(precision)
    params = np.asarray(rhos, dtype=dtype).reshape(-1, 1)
    y0 = np.zeros((len(params), 3), dtype=dtype)
    y0[:, 0] = 1
    problem = Lorenz()
    target = jax.devices(device)[0]
    atol = rtol / 1000

    if library == "GRADSOLVE":
        def run():
            with jax.default_device(target):
                result = gradsolve.solve(problem, y0, params, engine="cuda_tsit5",
                                         device=device, rtol=rtol, atol=atol)
            if result.solver != "cuda_tsit5":
                raise RuntimeError(f"GRADSOLVE rerouted away from CUDA Tsit5: {result.route}")
            return np.asarray(result.y_final)
    elif library == "Diffrax":
        def single(y, p):
            result = diffrax.diffeqsolve(
                diffrax.ODETerm(problem.f_jax), diffrax.Tsit5(),
                dtype.type(0), dtype.type(1), dtype.type(0.01), y, args=p,
                stepsize_controller=diffrax.PIDController(rtol=rtol, atol=atol),
                saveat=diffrax.SaveAt(t1=True), max_steps=100_000)
            return result.ys[0]
        solve = jax.jit(jax.vmap(single))

        def run():
            return np.asarray(solve(jax.device_put(y0, target),
                                    jax.device_put(params, target)))
    else:
        raise ValueError(f"Unknown library: {library}")

    def checked():
        result = run()
        if result.shape != y0.shape or result.dtype != dtype:
            raise AssertionError(f"Unexpected result shape/dtype: {result.shape}, {result.dtype}")
        if not np.isfinite(result).all():
            raise AssertionError("Non-finite solver output")
        return result

    return checked


def measure(run, samples=5):
    run()
    times = []
    for _ in range(samples):
        start = time.perf_counter()
        run()
        times.append(time.perf_counter() - start)
    return min(times)


def prepare_kernel(rhos, precision, rtol, device="gpu"):
    """Device-resident cuda_tsit5 kernel call and a host-side checked result.

    gradsolve exposes no public kernel-only API; this is the internal FFI runner
    its paper harness (benchmarks/forward_vs_diffeqgpu.py) times.
    """
    from gradsolve.cuda._ffi_bridge import make_runner

    dtype = np.dtype(precision)
    target = jax.devices(device)[0]
    y0 = np.zeros((3, len(rhos)), dtype=dtype)
    y0[0] = 1
    y0 = jax.device_put(y0, target)
    rho = jax.device_put(np.asarray(rhos, dtype=dtype), target)
    run = make_runner("lorenz", 3, precision, "cuda")

    def call():
        return jax.block_until_ready(run(y0, rho, 1.0, rtol, rtol / 1000, 100_000))

    def checked():
        result = np.ascontiguousarray(np.asarray(call()[0]).T)
        if result.shape != (len(rhos), 3) or result.dtype != dtype:
            raise AssertionError(f"Unexpected result shape/dtype: {result.shape}, {result.dtype}")
        if not np.isfinite(result).all():
            raise AssertionError("Non-finite solver output")
        return result

    return call, checked
