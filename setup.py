from setuptools import setup, find_packages

setup(
    name="pystran",
    version="0.1.0",  # Update with the actual version
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        # Add other dependencies here as needed
    ],
    author="Petr Krysl",
    author_email="pkrysl@ucsd.edu",
    description="Python STRuctural ANalysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/PetrKryslUCSD/pystran",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires=">=3.6",
)
