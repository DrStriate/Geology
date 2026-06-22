import pytest
from src.calculator import add, divide, subtract

def test_add_success():
    assert add(2, 3) == 5

def test_divide_success():
    assert divide(10, 2) == 5.0

def test_divide_by_zero_exception():
    with pytest.raises(ValueError, match="Cannot divide by zero."):
        divide(5, 0)

def test_subtract():
  assert subtract(2, 3) == -1
