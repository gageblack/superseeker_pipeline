from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='superseeker',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.19.0',
        'pandas>=1.1.0',
        'matplotlib>=3.3.0',
        'scipy>=1.5.0',
        'tqdm>=4.50.0',
    ],
    entry_points={
        'console_scripts': [
            'superseeker=superseeker.pipeline:run_pipeline',
        ],
    },
    author='Gage Black',
    author_email='gage.black@utah.edu',
    description='A Python library for identifying subclonal evolution in cancer',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/gageblack/superseeker_pipeline',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bioinformatics',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    include_package_data=True,
    license='MIT',
)
