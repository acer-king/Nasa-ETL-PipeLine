
import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="merra2",  # Replace with your own username
    version="0.0.1",
    author="Acer-king",
    author_email="acer-king03@gmail..com",
    description="First version Nasa ETL",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/acer-king/Nasa-ETL-PipeLine",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
