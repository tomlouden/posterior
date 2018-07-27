from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(	name='posterior', 
	version="0.1",
	author='Tom Louden',
	author_email = 'T.Louden@warwick.ac.uk',
	url = 'https://github.com/tomlouden/posterior',
	packages =['posterior'],
	license = ['GNU GPLv3'],
	description ='Just some simple tools for posterior inference',
	classifiers = [
		'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering',
		'Programming Language :: Python'
		],

)
