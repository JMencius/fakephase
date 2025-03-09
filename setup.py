from setuptools import setup, find_namespace_packages
setup(
    name = "fakephase",
    package_dir={"": "pie"},
    packages=find_namespace_packages(where="pie"),
    package_data={
        "pie.module": ["*.py"],
    },
    version = "0.1.0",
    description = "Fake telomere-to-telomere haloptype phasing algorithm",
    author = "Jun Mencius",
    author_email = "zjmeng22@m.fudan.edu.cn",
    url = "https://github.com/JMencius/fakephase",
    keywords = ["fakephase", "haplotype", "phasing"],
    python_requires = ">=3.10",
    install_requires = [
        "click>=8.1.8",
        "cyvcf2>=0.31.1",
        "pysam"
        ],
    extras_require = {
        "dev": ["pytest"],
        },
    entry_points={
    "console_scripts": [
        "pie = fakephase.fakephase:main",
        ],
    },
)
