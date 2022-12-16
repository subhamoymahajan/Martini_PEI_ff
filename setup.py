#    Copyright 2022 SUBHAMOY MAHAJAN 
#    
#    This file is part of InSilicoMicroscopy software.
# 
#    InSilicoMicroscopy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.)

import setuptools

setuptools.setup(
   name="Martini_PEI_ff",
   version="1.0.0 beta",
   author="Subhamoy Mahajan",
   author_email="subhamoy@ualberta.ca",
   description="Toolbox for generating and creating CG-PEI topology",
   url="https://github.com/subhamoymahajan/Martini_PEI_ff",
   license='GPLv3',
   install_requires=['networkx','numpy','scipy','matplotlib'],
   packages=['coarsen'],
   package_data={'': ['LICENSE.txt'],'coarsen': ['run_cg_sim.sh','cgff_2019.pickle','cgff_2022.pickle']},
   classifiers=[
        "Development Status :: 5-Production/Stable",
        "Intended Audience :: Science/Research",
        "Indended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language:: Python :: 3.7",
        "Programming Language:: Unix Shell",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
   ],
   python_requires='>=3.7',
   entry_points={
       'console_scripts': [
           'coarsen = coarsen.__main__:main'
       ]
   }
)
