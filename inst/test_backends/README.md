# `testthat` unit test suites for `MsBackend`

This folder contains unit test suites that can be called by other packages
defining `MsBackend` instances to ensure that they are compliant with the
`MsBackend` definition and that they provide all necessary
methods/functionality.

Different test suites are available within sub-folders of this directory, each
of them having one or multiple *.R* files defining the unit tests. A description
of the unit tests and which variables need to be defined by the calling package
is provided within these.

The unit tests from one test suite can be called from within another package (in
their `testthat.R` file) like shown in the example below that will call the unit
tests from the *test_MsBackend* test suite/directory:

```
test_suite <- system.file("test_backends", "test_MsBackend",
    package = "Spectra")
test_dir(test_suite)
```

## Authors/contributors

- Steffen Neumann (@sneumann)
- Helge Hecht (@hechth)
- Johannes Rainer (@jorainer)
