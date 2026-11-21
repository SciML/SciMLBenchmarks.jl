"""JAX Diffrax and PyTorch ensemble Lorenz timings.

Matches the ensemble Lorenz setup in Utkarsh et al., Comput. Methods Appl.
Mech. Eng. 428 (2024) 117109 (arXiv:2304.06835): Tsit5 / RK4, dt=0.001,
tspan=(0,1), rho in linspace(0, 21, N). Timing is done in Python so Julia
call overhead is not in the numbers. Min of repeats, after warmup, with a
device sync — the same statistic as the paper artifacts.
"""
import time

import jax
import jax.numpy as jnp
import numpy as np
import torch
import diffrax


def require_jax_gpu():
    devs = jax.devices("gpu")
    if not devs:
        raise RuntimeError("JAX did not see a GPU")
    return devs[0]


def require_torch_cuda():
    if not torch.cuda.is_available():
        raise RuntimeError("PyTorch did not see a CUDA GPU")
    return torch.device("cuda")


def _jax_array(values):
    return jax.device_put(np.asarray(values, dtype=np.float32), require_jax_gpu())


def _torch_array(values):
    return torch.tensor(values, dtype=torch.float32, device=require_torch_cuda())


def _lorenz(t, y, k1):
    f0 = 10.0 * (y[1] - y[0])
    f1 = k1 * y[0] - y[1] - y[0] * y[2]
    f2 = y[0] * y[1] - 2.666 * y[2]
    return jnp.stack([f0, f1, f2])


def _solve_diffrax(k1, adaptive):
    kwargs = {}
    if adaptive:
        kwargs["stepsize_controller"] = diffrax.PIDController(rtol=1e-8, atol=1e-8)
    sol = diffrax.diffeqsolve(
        diffrax.ODETerm(_lorenz),
        diffrax.Tsit5(),
        np.float32(0.0),
        np.float32(1.0),
        np.float32(0.001),
        jnp.array([1.0, 0.0, 0.0], dtype=jnp.float32),
        args=k1,
        saveat=diffrax.SaveAt(t0=True, t1=True),
        max_steps=100_000,
        **kwargs,
    )
    return sol.ys


_solve_fixed = jax.jit(jax.vmap(lambda k1: _solve_diffrax(k1, False)))
_solve_adaptive = jax.jit(jax.vmap(lambda k1: _solve_diffrax(k1, True)))


def _min_seconds(fn, warmup, repeats):
    for _ in range(warmup):
        fn()
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return float(min(times))


def time_diffrax(rhos, adaptive=False, warmup=3, repeats=5):
    rhos = _jax_array(rhos)
    solve = _solve_adaptive if adaptive else _solve_fixed

    def run():
        solve(rhos).block_until_ready()

    return _min_seconds(run, warmup, repeats)


def _lorenz_rhs_torch(u, rho):
    x = u[:, 0]
    y = u[:, 1]
    z = u[:, 2]
    du1 = 10.0 * (y - x)
    du2 = rho * x - y - x * z
    du3 = x * y - 2.666 * z
    return torch.stack((du1, du2, du3), dim=1)


def _rk4_batch(rho, dt=0.001, t1=1.0):
    # Same RK4 / dt as the paper's torchdiffeq method='rk4' run. Official
    # torchdiffeq does not vmap; the paper used a fork. Batched RK4 is the
    # same integrator on the same problem.
    steps = int(round(t1 / dt))
    u = torch.zeros(rho.shape[0], 3, device=rho.device, dtype=rho.dtype)
    u[:, 0] = 1.0
    u0 = u
    dt_t = torch.as_tensor(dt, device=rho.device, dtype=rho.dtype)
    for _ in range(steps):
        k1 = _lorenz_rhs_torch(u, rho)
        k2 = _lorenz_rhs_torch(u + 0.5 * dt_t * k1, rho)
        k3 = _lorenz_rhs_torch(u + 0.5 * dt_t * k2, rho)
        k4 = _lorenz_rhs_torch(u + dt_t * k3, rho)
        u = u + (dt_t / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return torch.stack((u0, u), dim=1)


def time_torch_rk4(rhos, warmup=2, repeats=5):
    rho = _torch_array(rhos)

    def run():
        _rk4_batch(rho)
        torch.cuda.synchronize()

    return _min_seconds(run, warmup, repeats)


def _checked_host(value, shape):
    if isinstance(value, torch.Tensor):
        if not value.is_cuda:
            raise AssertionError("PyTorch output is not on CUDA")
        host = value.detach().cpu().numpy()
    else:
        if not all(device.platform == "gpu" for device in value.devices()):
            raise AssertionError("JAX output is not on a GPU")
        host = np.asarray(value)
    if host.shape != shape or host.dtype != np.float32:
        raise AssertionError(f"Unexpected shape/dtype: {host.shape}, {host.dtype}")
    if not np.isfinite(host).all():
        raise AssertionError("Non-finite output")
    return host


def checked_outputs(rhos, state, rho):
    """Return untimed GPU input, RHS, and solution checks as host Float32 arrays."""
    jax_rhos, torch_rhos = _jax_array(rhos), _torch_array(rhos)
    jax_rhs = _lorenz(np.float32(0), _jax_array(state), np.float32(rho))
    torch_rhs = _lorenz_rhs_torch(_torch_array(state)[None, :], _torch_array([rho]))[0]
    n = len(rhos)
    return {
        "JAX parameters": _checked_host(jax_rhos, (n,)),
        "PyTorch parameters": _checked_host(torch_rhos, (n,)),
        "JAX RHS": _checked_host(jax_rhs, (3,)),
        "PyTorch RHS": _checked_host(torch_rhs, (3,)),
        "Diffrax fixed": _checked_host(_solve_fixed(jax_rhos), (n, 2, 3)),
        "Diffrax adaptive": _checked_host(_solve_adaptive(jax_rhos), (n, 2, 3)),
        "PyTorch RK4": _checked_host(_rk4_batch(torch_rhos), (n, 2, 3)),
    }
