
"""
statreporter.py:  Provides capabilities to report GaMD statistics

Portions copyright (c) 2021 University of Kansas
Authors: Matthew Copeland

"""


class StatisticsReporter:

    def __init__(self, recording_rate, filename):
        self.recording_rate = recording_rate
        self.filename = filename
        self.lastValues = []

    def report(self, step, integrator):
        pass
