class Klass:
    def __init__(self, i):
        global z
        z = i+1
        self.i = i

    def p(self):
        print(z)
        
