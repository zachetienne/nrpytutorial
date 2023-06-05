# https://stackoverflow.com/questions/3393612/run-certain-code-every-n-seconds/13151299

from threading import Timer

class RepeatedTimer(object):
    """
    A Timer class that can executes a function every 'interval' seconds.
    """

    def __init__(self, interval, function, *args, **kwargs):
        """
        Creates a timer.

        :param interval: int or float, the number of seconds between each call to 'function'
        :param function: function, the function to call
        :param args: list, optional, positional arguments to pass to 'function'
        :param kwargs: dict, optional, keyword arguments to pass to 'function'
        """
        self._timer     = None
        self.interval   = interval
        self.function   = function
        self.args       = args
        self.kwargs     = kwargs
        self.is_running = False
        self.start()

    def _run(self):
        """
        Is called every 'interval' seconds and runs 'function'.
        """
        self.is_running = False
        self.start()
        try:
            self.function(*self.args, **self.kwargs)
        except Exception as e:
            # You could add better exception handling here
            print(str(e))

    def start(self):
        """
        Starts the timer.
        """
        if not self.is_running:
            self._timer = Timer(self.interval, self._run)
            self._timer.start()
            self.is_running = True

    def stop(self):
        """
        Stops the timer.
        """
        self._timer.cancel()
        self.is_running = False
