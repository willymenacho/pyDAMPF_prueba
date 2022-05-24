import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyDAMPF",
    version="0.0.1",
    author="Horacio V. Guzman & Willy Menacho N.",
    author_email="horacio.guzman@ijs.si",
    description="pyDAMPF is a tool oriented to the Atomic Force Microscopy (AFM) community, which allows the simulation of the physical properties of materials under variable relative humidity (RH).",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/willymenacho/pyDAMPF_prueba",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)
