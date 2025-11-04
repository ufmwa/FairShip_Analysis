import contextlib

try:
    import torch.nn as nn
    _TORCH_AVAILABLE = True
except Exception:
    nn = None
    _TORCH_AVAILABLE = False


if _TORCH_AVAILABLE:

    class AbstractModule(nn.Module):
        def __init__(self, *args, **kwargs):
            super(AbstractModule, self).__init__()

        @contextlib.contextmanager
        def _enter_variable_scope(self, *args, **kwargs):
            yield None


else:

    class AbstractModule:
        def __init__(self, *args, **kwargs):
            raise RuntimeError(
                "AbstractModule requires torch; caller should have gated this path."
            )

        @contextlib.contextmanager
        def _enter_variable_scope(self, *args, **kwargs):
            yield None
