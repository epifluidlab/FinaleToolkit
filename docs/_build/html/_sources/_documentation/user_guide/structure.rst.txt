
Structure
=========================================

While not particularly pertinent to your usage, we include this page to give you a sense of what is going on under the hood of **FinaleToolkit**:

At the heart of **FinaleToolkit** is the ``Fragment Generator``. This generator is designed to be as flexible as possible. It takes BAM, CRAM, SAM, and Fragment Files (among other arguments) and filters out the reads that are not of interest (e.x. unmapped reads, reads with low mapping quality, etc.). The fragments that make it through this filter are then passed to the relevant function to be processed.

A rough outline of the structure of **FinaleToolkit** is therefore as follows:

.. figure:: finaletoolkit_structure.png
   :align: center