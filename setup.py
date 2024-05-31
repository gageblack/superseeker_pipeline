from setuptools import setup, find_packages

setup(
    name='superseeker',
    version='0.1.dev3',
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
    author='Gage Black',
    author_email='gage.black@utah.edu',
    description='SuperSeeker_Pipeline is a Python library for identifying subclonal evolution in cancer.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/gageblack/superseeker_pipeline',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    include_package_data=True,
    license='MIT',
)
