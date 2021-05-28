using Gtk

function f(cos)
    show(time)
end

function d(cos)
    show(freq_window)
end

b = GtkBuilder(filename="gui1.glade")
win = b["window"]
time = b["time_window"]
freq_window = b["freq_window"]
btn_time = b["btn_time"]
btn_freq = b["btn_freq"]

signal_connect(f,btn_time,"clicked")
signal_connect(d,btn_freq,"clicked")


showall(win)