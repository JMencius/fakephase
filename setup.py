from setuptools import setup, find_namespace_packages
setup(
    name = "fakephase",
    package_dir={"": "src"},
    packages=find_namespace_packages(where="src"),
    package_data={
        "fakephase.modules": ["*.py"],
        "fakephase.classes": ["*.py"],
    },
    version = "0.1.0",
    description = "Fake telomere-to-telomere haloptype phasing algorithm",
    author = "Jun Mencius",
    author_email = "zjmeng22@m.fudan.edu.cn",
    url = "https://github.com/JMencius/fakephase",
    keywords = ["fakephase", "haplotype", "phasing"],
    python_requires = ">=3.7",
    install_requires = [
        "click>=8.1.8",
        "cyvcf2>=0.31.1",
        "pyfastx>=2.2.0",
        "pysam>=0.23.0",
        ],
    extras_require = {
        "dev": ["pytest"],
        },
    entry_points={
    "console_scripts": [
        "fakephase = fakephase.fakephase:main",
        ],
    },
)
