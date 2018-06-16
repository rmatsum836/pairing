"""
pairing
analyze pairing and clustering of molecular systems
"""
from setuptools import setup
import versioneer

DOCLINES = __doc__.split("\n")

setup(
    name='pairing',
    author='Matthew W. Thompson',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',
    packages=['pairing', "pairing.tests"],
    platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
    python_requires=">=3.5",
    zip_safe=False,
)
