# pyDAMPF_prueba
prueba para diversos fines 

pyDAMPF: a Python package for modeling mechanical properties of hygroscopic materials under interaction with a nanoprobe
======================================================

|CI Status| |Coverage Status| |Documentation Status|

[//]: # ( .. |CI Status| image:: https://github.com/pypr/compyle/actions/workflows/tests.yml/badge.svg)
[//]: # (    :target: https://github.com/pypr/compyle/actions/workflows/tests.yml)
[//]: # ( .. |Documentation Status| image:: https://readthedocs.org/projects/compyle/badge/?version=latest)
[//]: # (  :target: https://compyle.readthedocs.io/en/latest/?badge=latest)
[//]: # ( :alt: Documentation Status)
[//]: # ( .. |Coverage Status| image:: https://codecov.io/gh/pypr/compyle/branch/master/graph/badge.svg)
[//]: # ( :target: https://codecov.io/gh/pypr/compyle)

pyDAMPF is a tool oriented to the Atomic Force Microscopy (AFM) community, which allows the simulation of the physical properties of materials under variable relative humidity (RH).

The computing engine is written in Fortran in order to reuse physics code and wrapped to Python to interoperate with high-level packages. We also introduce an in-house multi-thread approach of the pyDAMPF code, which compares for various computing architectures (PC, Google Colab and a HPC facility) very favorable in comparison to a former AFM simulator. 


Users can use the existing cantilever database to perform their simulations or add-up further models. There are three ways of running pyDAMPF:

- serial (uses a single thread)
- multi-thread (unix based computer)
- multi-thread (SLURM based computer)


Documentation and learning material is also available in the form of:

- An introduction to compyle in the context of writing a parallel molecular
  dynamics simulator is in our `SciPy 2020 paper
  <http://conference.scipy.org/proceedings/scipy2020/compyle_pr_ab.html>`_.

- `Compyle poster presentation <https://docs.google.com/presentation/d/1LS9XO5pQXz8G5d27RP5oWLFxUA-Fr5OvfVUGsgg86TQ/edit#slide=id.p>`_

Installation
-------------

Compyle is itself largely pure Python but depends on numpy_ and requires
either Cython_ or PyOpenCL_ or PyCUDA_ along with the respective backends of a
C/C++ compiler, OpenCL and CUDA. If you are only going to execute code on a
CPU then all you need is Cython.

You should be able to install Compyle by doing::

  $ pip install compyle
