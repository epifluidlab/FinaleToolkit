==========================
Intersect Policy
==========================

When you give the ``Fragment Generator`` an interval of ``start`` and ``end`` to pull fragments from, which fragments will it return? 

What if a fragment has 1bp of overlap within the just on the edge of the window? Does it also return that fragment as well?

These questions are answered by the ``intersect_policy`` argument, which shows up in many of the CLI and API commands.

``intersect_policy`` can either be ``midpoint`` or ``any``.

Midpoint
----------------

``midpoint``: The midpoint of the fragment determines whether the fragment will be returned or not. If the midpoint falls between ``start`` and ``end``, then it will be returned. Else, it will not.

Any
-------

``any``: As long as the fragment overlaps the interval between ``start`` and ``end`` in some capacity, then it will be returned. Else, it will not.