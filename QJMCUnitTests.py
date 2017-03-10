import unittest
import QJMCMathUnitTests
import nose

print('MATH TESTS')
nose.run(argv=[__file__, 'QJMCMathUnitTests.py'])
print('SETUP TESTS')
nose.run(argv=[__file__, 'QJMCSetUpUnitTests.py'])
