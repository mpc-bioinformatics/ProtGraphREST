class ProtGraphException(Exception):
    def __init__(self, status, body, headers={"content-type": "application/json"}):
        self.status = status
        self.body = body
        self.headers = headers
