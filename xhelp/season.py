class Season:
    def __init__(self, DISCRETE_GENERATIONS):
        self.countdown = float("inf")
        self.DISCRETE_GENERATIONS = DISCRETE_GENERATIONS
        self.reset()

    def reset(self):
        self.countdown = (
            self.DISCRETE_GENERATIONS if self.DISCRETE_GENERATIONS else float("inf")
        )
