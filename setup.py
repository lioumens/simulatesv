from setuptools import setup, find_packages

setup(
        name='simulatesv',
        version='0.2.0',
        description='Simulate structural variations and SNPs with artificial dna sequences',
        long_description='''Simulatesv is a package to generate artificial
        genomes and some variants. Simulatesv is used to generate a validation
        set and quick testing tool when evaluating variant detection tools.''',
        url='https://github.com/mlliou112/simulatesv',
        author='Michael Liou',
        author_email='mliou112@yahoo.com',
        license='MIT',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Topic :: Scientific/Engineering :: Bio-Informatics'
            ],
        keywords='bioinformatics structural-variation simulate dna genetics sequencing',
        packages=find_packages(),
        install_requires=['numpy'],
        test_suite='nose.collector',
        tests_require=['nose'],
        scripts=['simulatesv/simulatesv.py'])
