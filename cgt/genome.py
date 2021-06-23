"""
"""

class Genome:
    """Optional wrapper for genomes, so that they can be printed/operated on in a more pythonic way."""
    def __init__(self, framework, genome):
        self.framework = framework
        self.generating_object = genome
        try:
            self.canonical_instance = framework.canonical_instance(framework(genome))
        except:
            raise TypeError("Unable to create genome object with given data")