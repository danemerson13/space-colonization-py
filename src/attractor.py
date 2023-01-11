class Attractor:
    # Properties of an attractor:
    # Location: type = 3x1 double array
    # Target: type = Boolean, if target attractor, cannot be killed until reached
    def __init__(self,loc,target):
        self.location = loc
        self.target = target
        self.killed = False

    def getLocation(self):
        return self.location

    def isTarget(self):
        return self.target

    def setKill(self):
        self.killed = True

    def isKilled(self):
        return self.killed