"""
Python timing utilities and model definitions for neural network benchmarks.
All timing is done inside Python to avoid Julia-to-Python call overhead.
"""
import time
import numpy as np

# ---- JAX ----
import jax
import jax.numpy as jnp
from jax import random
import optax

# ---- PyTorch ----
import torch
import torch.nn as tnn


# ============================================================
# Timing utilities
# ============================================================

def time_jax_inference(fn, *args, n_runs=100):
    """Time a JIT-compiled JAX function, including block_until_ready."""
    for _ in range(5):
        fn(*args).block_until_ready()
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        fn(*args).block_until_ready()
        end_ = time.perf_counter()
        times.append(end_ - start)
    return float(np.median(times))


def time_jax_train_step(train_step_fn, params, x, y, n_runs=50):
    """Time a JAX training step, blocking until all outputs are ready."""
    for _ in range(3):
        params, _ = train_step_fn(params, x, y)
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        params, loss = train_step_fn(params, x, y)
        jax.tree.map(lambda t: t.block_until_ready(), params)
        end_ = time.perf_counter()
        times.append(end_ - start)
    return float(np.median(times))


def time_torch_inference(model, x, n_runs=100):
    """Time a PyTorch model in eval mode."""
    model.eval()
    with torch.no_grad():
        for _ in range(5):
            model(x)
        times = []
        for _ in range(n_runs):
            start = time.perf_counter()
            model(x)
            end_ = time.perf_counter()
            times.append(end_ - start)
    return float(np.median(times))


def time_torch_train_step(model, optimizer, loss_fn, x, y, n_runs=50):
    """Time a full PyTorch train step (zero_grad + forward + backward + step)."""
    model.train()
    for _ in range(3):
        optimizer.zero_grad()
        out = model(x)
        loss = loss_fn(out, y)
        loss.backward()
        optimizer.step()
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        optimizer.zero_grad()
        out = model(x)
        loss = loss_fn(out, y)
        loss.backward()
        optimizer.step()
        end_ = time.perf_counter()
        times.append(end_ - start)
    return float(np.median(times))


# ============================================================
# JAX MLP helpers
# ============================================================

def create_jax_mlp_params(key):
    """Create parameters for a 7-layer MLP (32 -> 256x6 -> 10)."""
    params = []
    key, subkey = random.split(key)
    params.append({"w": random.normal(subkey, (32, 256)) * 0.01, "b": jnp.zeros(256)})
    for _ in range(5):
        key, subkey = random.split(key)
        params.append(
            {"w": random.normal(subkey, (256, 256)) * 0.01, "b": jnp.zeros(256)}
        )
    key, subkey = random.split(key)
    params.append({"w": random.normal(subkey, (256, 10)) * 0.01, "b": jnp.zeros(10)})
    return params


@jax.jit
def jax_mlp_relu_forward(params, x):
    for layer in params[:-1]:
        x = jax.nn.relu(x @ layer["w"] + layer["b"])
    return x @ params[-1]["w"] + params[-1]["b"]


@jax.jit
def jax_mlp_gelu_forward(params, x):
    for layer in params[:-1]:
        x = jax.nn.gelu(x @ layer["w"] + layer["b"])
    return x @ params[-1]["w"] + params[-1]["b"]


def make_jax_mlp_train_step(forward_fn, params):
    """Create a JIT-compiled training step for a JAX MLP."""
    optimizer = optax.adam(1e-3)
    opt_state = optimizer.init(params)

    def loss_fn(p, x, y):
        pred = forward_fn(p, x)
        return jnp.mean((pred - y) ** 2)

    @jax.jit
    def train_step(params, x, y):
        loss, grads = jax.value_and_grad(loss_fn)(params, x, y)
        updates, _ = optimizer.update(grads, opt_state, params)
        new_params = optax.apply_updates(params, updates)
        return new_params, loss

    return train_step


# ============================================================
# JAX LeNet helpers
# ============================================================

def create_jax_lenet_params(key):
    params = {}
    key, k1 = random.split(key)
    params["conv1_w"] = random.normal(k1, (6, 1, 5, 5)) * 0.01
    params["conv1_b"] = jnp.zeros(6)
    key, k2 = random.split(key)
    params["conv2_w"] = random.normal(k2, (16, 6, 5, 5)) * 0.01
    params["conv2_b"] = jnp.zeros(16)
    key, k3 = random.split(key)
    params["fc1_w"] = random.normal(k3, (256, 120)) * 0.01
    params["fc1_b"] = jnp.zeros(120)
    key, k4 = random.split(key)
    params["fc2_w"] = random.normal(k4, (120, 84)) * 0.01
    params["fc2_b"] = jnp.zeros(84)
    key, k5 = random.split(key)
    params["fc3_w"] = random.normal(k5, (84, 10)) * 0.01
    params["fc3_b"] = jnp.zeros(10)
    return params


@jax.jit
def jax_lenet_forward(params, x):
    """LeNet-5 forward pass. x: (batch, 1, 28, 28) NCHW."""
    x = (
        jax.lax.conv(x, params["conv1_w"], (1, 1), "VALID")
        + params["conv1_b"][None, :, None, None]
    )
    x = jax.nn.relu(x)
    x = jax.lax.reduce_window(
        x, -jnp.inf, jax.lax.max, (1, 1, 2, 2), (1, 1, 2, 2), "VALID"
    )
    x = (
        jax.lax.conv(x, params["conv2_w"], (1, 1), "VALID")
        + params["conv2_b"][None, :, None, None]
    )
    x = jax.nn.relu(x)
    x = jax.lax.reduce_window(
        x, -jnp.inf, jax.lax.max, (1, 1, 2, 2), (1, 1, 2, 2), "VALID"
    )
    x = x.reshape(x.shape[0], -1)
    x = jax.nn.relu(x @ params["fc1_w"] + params["fc1_b"])
    x = jax.nn.relu(x @ params["fc2_w"] + params["fc2_b"])
    return x @ params["fc3_w"] + params["fc3_b"]


def make_jax_lenet_train_step(params):
    optimizer = optax.adam(1e-3)
    opt_state = optimizer.init(params)

    def loss_fn(p, x, y):
        pred = jax_lenet_forward(p, x)
        return jnp.mean((pred - y) ** 2)

    @jax.jit
    def train_step(params, x, y):
        loss, grads = jax.value_and_grad(loss_fn)(params, x, y)
        updates, _ = optimizer.update(grads, opt_state, params)
        return optax.apply_updates(params, updates), loss

    return train_step


# ============================================================
# PyTorch model definitions
# ============================================================

class TorchMLP(tnn.Module):
    def __init__(self, activation="relu"):
        super().__init__()
        act_fn = tnn.ReLU() if activation == "relu" else tnn.GELU()
        layers = [tnn.Linear(32, 256), act_fn]
        for _ in range(5):
            layers.extend([tnn.Linear(256, 256), act_fn])
        layers.append(tnn.Linear(256, 10))
        self.net = tnn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


class TorchMLPBN(tnn.Module):
    def __init__(self):
        super().__init__()
        layers = [tnn.Linear(32, 256), tnn.BatchNorm1d(256), tnn.ReLU()]
        for _ in range(5):
            layers.extend([tnn.Linear(256, 256), tnn.BatchNorm1d(256), tnn.ReLU()])
        layers.append(tnn.Linear(256, 10))
        self.net = tnn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


class TorchLeNet(tnn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = tnn.Conv2d(1, 6, 5)
        self.conv2 = tnn.Conv2d(6, 16, 5)
        self.pool = tnn.MaxPool2d(2, 2)
        self.fc1 = tnn.Linear(16 * 4 * 4, 120)
        self.fc2 = tnn.Linear(120, 84)
        self.fc3 = tnn.Linear(84, 10)

    def forward(self, x):
        x = self.pool(torch.relu(self.conv1(x)))
        x = self.pool(torch.relu(self.conv2(x)))
        x = x.view(x.size(0), -1)
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        return self.fc3(x)


class TorchResBlock(tnn.Module):
    def __init__(self, channels):
        super().__init__()
        self.conv1 = tnn.Conv2d(channels, channels, 3, padding=1, bias=False)
        self.bn1 = tnn.BatchNorm2d(channels)
        self.conv2 = tnn.Conv2d(channels, channels, 3, padding=1, bias=False)
        self.bn2 = tnn.BatchNorm2d(channels)

    def forward(self, x):
        out = torch.relu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        return torch.relu(out + x)


class TorchSmallResNet(tnn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = tnn.Conv2d(3, 64, 3, padding=1, bias=False)
        self.bn1 = tnn.BatchNorm2d(64)
        self.block1 = TorchResBlock(64)
        self.down_conv1 = tnn.Conv2d(64, 128, 3, stride=2, padding=1, bias=False)
        self.down_bn1 = tnn.BatchNorm2d(128)
        self.down_conv2 = tnn.Conv2d(128, 128, 3, padding=1, bias=False)
        self.down_bn2 = tnn.BatchNorm2d(128)
        self.pool = tnn.AdaptiveAvgPool2d(1)
        self.fc = tnn.Linear(128, 10)

    def forward(self, x):
        x = torch.relu(self.bn1(self.conv1(x)))
        x = self.block1(x)
        x = torch.relu(self.down_bn1(self.down_conv1(x)))
        x = self.down_bn2(self.down_conv2(x))
        x = self.pool(x).view(x.size(0), -1)
        return self.fc(x)
