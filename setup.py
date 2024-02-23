import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='comparative_genomics',
    version='0.17',
    packages=setuptools.find_packages(where='src'),
    url='https://github.com/kinestetika/comparative_genomics',
    license='MIT',
    author='Marc Strous',
    author_email='mstrous@ucalgary.ca',
    description='Call orthologues and create a phylogenetic tree based on concatenated alignment of conserved proteins',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Natural Language :: English',
                 'Operating System :: OS Independent',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3.10',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='taxonomy phylogeny alignment',
    project_urls={'Source': 'https://github.com/kinestetika/comparative_genomics'},
    package_dir={'': 'src'},
    package_data={'': ['gtdb-pfam.hmm', 'gtdb-tigr.hmm', 'ribosomal.pfam.hmm', 'ribosomal.tigr.hmm', 'rpoABC.hmm']},
    python_requires='>=3.10',
    install_requires=['clipkit'],
    extras_require={  # Optional
        'dev': ['setuptools', 'build', 'twine'],
        'test': []
    },
    entry_points={  # Optional
        'console_scripts': [
            'tree_of_mags=comparative_genomics.tree_of_mags:main',
            'orthologues=comparative_genomics.orthologues:main',
        ],
    }
)