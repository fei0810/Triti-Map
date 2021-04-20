from os.path import dirname, join, abspath
import setuptools

this_directory = abspath(dirname(__file__))
with open(join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setuptools.setup(
    name="tritimap",
    version="0.9.4",
    author="Fei Zhao",
    author_email="zhaofei920810@gmail.com",
    url="https://github.com/fei0810/Triti-Map",
    description="A Snakemake-based pipeline for gene mapping in Triticeae.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['tritimap'],
    package_data={'': ["tritimap/*", ]},
    data_files=[(".", ["README.md"])],
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        'click',
        'snakemake'
    ],
    entry_points={
        'console_scripts': [
            'tritimap = tritimap.tritimap:cli'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Intended Audience :: Science/Research',
        'Development Status :: 4 - Beta',
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
