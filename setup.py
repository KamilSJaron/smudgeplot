import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="smudgeplot-KamilSJaron",
    version="0.1.13.dev0",
    author="Kamil S. Jaron, Timothy Rhyker Ranallo-Benavidez",
    author_email="kamiljaron+smudge@gmail.com",
    description="a package for inference of ploidy and heterozygosity structure using whole genome sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tbenavi1/smudgeplot",
    packages=setuptools.find_packages(),
    scripts=['scripts/smudgeplot'],
    zip_safe=False,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
