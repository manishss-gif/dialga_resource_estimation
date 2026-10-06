import unittest
import test_ct_primitives
import test_ct_profiles

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromModule(test_ct_primitives))
    suite.addTests(loader.loadTestsFromModule(test_ct_profiles))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
