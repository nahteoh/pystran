from setuptools import setup, find_packages

setup(
    name="pystran",
    version="0.1.0",  # Replace with the actual version of your package
    author="Petr Krysl",
    author_email="pkrysl@ucsd.edu",  # Replace with the author's actual email if provided
    description="A Python package for structural analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/PetrKryslUCSD/pystran",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Adjust if the license is different
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    install_requires=[
        # Include package dependencies here
        "numpy",
        "scipy",
        "matplotlib",  # Add or remove based on the project's requirements
    ],
    entry_points={
    },
)
