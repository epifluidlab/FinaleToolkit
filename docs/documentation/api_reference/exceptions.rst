Exceptions
======================================

FinaleToolkit raises informative, typed exceptions. Each one also subclasses the
built-in exception it replaces (``ValueError``, ``FileNotFoundError``, or
``IndexError``), so existing ``except ValueError:`` / ``except FileNotFoundError:``
handlers continue to catch them.

.. automodule:: finaletoolkit.exceptions
   :members:
   :show-inheritance:
