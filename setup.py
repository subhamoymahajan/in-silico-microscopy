#    Copyright 2020,2021 SUBHAMOY MAHAJAN 
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
   name="in-silico-microscopy",
   version="1.2.1",
   author="Subhamoy Mahajan",
   author_email="subhamoygithub@gmail.com",
   description="Toolbox for generating in-silico microscopy images from molecular simulations",
   url="https://github.com/subhamoymahajan/in-silico-microscopy",
   license='GPLv3',
   install_requires=['numpy','matplotlib','opencv-python'],
   packages=['siliscopy'],
   package_data={'': ['LICENSE.txt'],'siliscopy': ['gen_mono', 'gen_mono.c']},
   classifiers=[
        "Development Status :: 5-Production/Stable",
        "Intended Audience :: Science/Research",
        "Indended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language:: Python :: 3.7",
        "Programming Language:: C",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
   ],
   python_requires='>=3.7',
   entry_points={
       'console_scripts': [
           'siliscopy = siliscopy.__main__:main'
       ]
   }
)
