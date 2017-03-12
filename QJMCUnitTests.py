#Edited 12/3/17 Ben Everest
#A collection of all the unit tests, just run this to run all
import unittest
import QJMCMathUnitTests
import nose

print('MATH TESTS')
nose.run(argv=[__file__, 'QJMCMathUnitTests.py'])
print('SETUP TESTS')
nose.run(argv=[__file__, 'QJMCSetUpUnitTests.py'])
