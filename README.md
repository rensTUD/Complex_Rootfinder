# Complex Rootfinder

A Python package for finding complex roots using the argument principle.

## Installation

```bash
pip install -e .
```

## Usage

```python
from complex_rootfinder.argument_principle import argument_principle

# Example usage
def f(z):
    return z**4 - 16  # finds roots at ±2 and ±2i

roots = argument_principle(
    real_min=-4,
    real_max=4,
    imag_min=-4,
    imag_max=4,
    step_size=0.1,
    det_func=f
)
```

Now update the import in your test file:

```python:tests/tests.py
from typing import Callable, Tuple
import numpy as np
import sys
import os

from complex_rootfinder.argument_principle import argument_principle

# ... rest of the test file remains the same ...
```

Also, create an empty `__init__.py` file in your package directory:

```python:complex_rootfinder/__init__.py
"""Complex rootfinder package."""
```

Your project structure should look like this: 
