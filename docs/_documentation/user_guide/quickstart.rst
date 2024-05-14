
Quickstart
=========================================

------------------------
Getting Started    
------------------------

If not to be integrated into a workflow, **FinaleToolkit** is intended to be run directly from a terminal interface. 

You can open the terminal and typing the following command::

    $ finaletoolkit

------------------------
Importing Modules   
------------------------

The **FinaleToolkit** package itself is divided by its API functions.

You can load *all* of the functions by using the following command::

    >>> import finaletoolkit as ft
    
If you want to load a *specific* function, you can do so by using the following:

Assuming we wanted to load the DELFI function::

    >>> from finaletoolkit.frag import delfi
