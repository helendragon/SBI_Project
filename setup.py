from setuptools import setup, find_packages 


setup(name='complexbuilder', 
    version='1.0', 
    description='Build a macrocomplex from interaction pairs',
	author='Nuria Olvera and Helena Catena',
	author_email='nuriaolvera17@gmail.com, helena.catena01@estudiant.upf.edu',
	url='https://github.com/helendragon/SBI_Project.git',
    packages= ['folder' ], 
    scripts=['scripts/complexbuilder'],
    install_requires=['biopython' ]
    )