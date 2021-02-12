import setuptools

setuptools.setup(
   name="in-silico-microscopy",
   version="1.2.0",
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
