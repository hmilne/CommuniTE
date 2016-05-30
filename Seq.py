class Seq:
    def __init__(self, num, seqID):
        self.num = num
        self.ID = seqID
        # self.sequence = sequence
        self.attributes = {}
    def setAttribute(self, key, value):
        self.attributes[key] = value
    def length(self):
        return len(self.sequence)
    def setSequence(self, sequence):
        self.sequence = sequence




