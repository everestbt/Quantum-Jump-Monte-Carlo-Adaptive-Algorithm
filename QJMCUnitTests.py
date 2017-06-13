#Edited 12/3/17 Ben Everest
#A collection of all the unit tests, just run this to run all
import nose

print('MATH TESTS')
nose.run(argv=[__file__, 'QJMCMathUnitTests.py'])
print('SETUP TESTS')
nose.run(argv=[__file__, 'QJMCSetUpUnitTests.py'])
print('MEASURE TESTS')
nose.run(argv=[__file__, 'QJMCMeasureUnitTests.py'])
print('JUMP TESTS')
nose.run(argv=[__file__, 'QJMCJumpUnitTests.py'])
