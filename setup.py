from setuptools import setup, find_packages

setup(
    name='SuperSeeker',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
    ],
    entry_points={
        'console_scripts': [
            'superseeker=superseeker.pipeline:run_pipeline',
        ],
    },
    author='Your Name',
    author_email='your.email@example.com',
    description='SuperSeeker is a Python library for identifying subclonal evolution in cancer.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/gageblack/superseeker_pipeline',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    include_package_data=True,
    license='MIT',
)