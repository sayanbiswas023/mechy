from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.2'
DESCRIPTION = 'Core Engineering Analysis Functionalities'

# Setting up
setup(
    name="mechy",
    version=VERSION,
    author="Sayan Biswas",
    author_email="<biswas0235sayan@proton.me>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=['numpy','matplotlib'],
    keywords=['python', 'engineering','finite element analysis','mechanics','stress','strain'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
