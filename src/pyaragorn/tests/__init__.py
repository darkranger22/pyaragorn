from . import (
    test_rna_finder,
    test_doctest,
)


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_rna_finder))
    test_doctest.load_tests(loader, suite, pattern)
    return suite
