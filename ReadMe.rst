WATER_BUDGET
============

``@author`` Marialaura Bancheri

``@copyright`` GNU Public License v3 AboutHydrology (Riccardo Rigon)

Description
-----------

This is a Gradle project

.. figure:: doc/ReadMe/gradle.png
   :alt: Gradle logo

   Gradle logo

To build the project just run

::

   ~ $ gradle build

in the working directory. You will find the built ``.jar`` in
``build/libs``. But a built version of the lates release is already
available at the `GitHub release
section <release://github.com/geoframecomponents/WATER_BUDGET/releases>`__

To build the ReadMe file from the markdown one in the doc/ReadMe folder

::

   pandoc doc/ReadMe/ReadMe.md -o ReadMe.rst --bibliography=doc/ReadMe/biblio.bib

Implementation
~~~~~~~~~~~~~~

Documentation
-------------

`Model Equation <doc/ReadMe/Model_equation.md>`__

Canopy
------

Canopy interception model (Rutter 1971; Aston 1979)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Developers’ documentation
~~~~~~~~~~~~~~~~~~~~~~~~~

Linkers’ documentation
~~~~~~~~~~~~~~~~~~~~~~

The documentation is available at
http://www.slideshare.net/marialaurabancheri/jgrassnewage-water-budget

Users’ documentation
~~~~~~~~~~~~~~~~~~~~

References
----------
