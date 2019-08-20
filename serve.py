# create the server
# noinspection PyCompatibility
from http.server import HTTPServer, SimpleHTTPRequestHandler
port = 8000
def serve(server_class=HTTPServer, handler_class=SimpleHTTPRequestHandler):
    server_address = ('', port)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()

# start the server in a separate thread to avoid blocking the main thread, thereby keeping it responsive.
# daemon threads in python: https://realpython.com/intro-to-python-threading/
import threading
threading.Thread(target=serve, args=(), daemon=True).start()

# open the browser
import webbrowser
webbrowser.open_new('http://localhost:{}'.format(port))

# stop server when ctrl-C is pressed
import sys
import time
while True:
    try:
        time.sleep(1)
    except KeyboardInterrupt:
        sys.exit(0)
