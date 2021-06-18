using Plots: length
using Gtk
using WAV
using Plots
using FFTW

include("naszmodul.jl")#ściezka do pliku naszmodul.jl

file = Nothing

sound = Nothing
fs = Nothing
sound_fftw = Nothing
sound_freq = Nothing
time_range = Nothing

temporary_sound = Nothing
temporary_fs = Nothing
temporary_fftw = Nothing
temporary_time_range =Nothing
temporary_freq = Nothing

current_sound = Nothing
current_fs = Nothing
current_fftw = Nothing
current_time_range =Nothing
current_freq = Nothing
current2_sound = Nothing
current2_fftw = Nothing
current3_sound = Nothing
#wczytywanie pliku
function file_load(cos)
    global file = open_dialog("Wybierz plik",GtkNullContainer(),("*.wav",))
    sound, fsp = wavread(file)

    global sound = sound[:, 1]
    global temporary_sound = sound
    global fs = fsp
    global temporary_fs = fs
    
    global sound_fftw = fft(sound)
    global temporary_fftw = sound_fftw
    global sound_freq = fftfreq(length(sound),fs)
    global temporary_freq = sound_freq
    global time_range = LinRange(0,length(sound)/fs,length(sound))
    global temporary_time_range = time_range

    p = Plots.plot(time_range,sound,xlabel="czas",ylabel = "amplituda", label= "nagranie")
    Plots.savefig(p,"plot_time.png")

    set_gtk_property!(main_plot,:file,"plot_time.png")
    set_gtk_property!(plot_time_menu,:file,"plot_time.png")
    set_gtk_property!(adjustment4,:upper,length(sound)/fs)
    set_gtk_property!(adjustment3,:upper,length(sound)/fs)
    set_gtk_property!(adjustment4,:value,length(sound)/fs)
end
#okno dotyczące czasu
function time_menu_open(cos) 
    set_gtk_property!(main_plot,:file,"plot_time.png")
    global current_sound = temporary_sound
    global current2_sound = temporary_sound
    global current_fs = temporary_fs
    global current_time_range = temporary_time_range
    show(time)
end
function open_time_cut(cos)
    set_gtk_property!(time_cut_image,:file,"plot_time.png")
    show(time_cut_window)
end
function cut_time(cos)
    start = get_gtk_property(adjustment3, :value, Float64)
    stop = get_gtk_property(adjustment4,:value,Float64)
    if start < stop
        time_new, sound_new = NaszModul.cutting_time(current_sound,current_fs,start,stop)
        global current2_sound = sound_new
        global current2_time_range = time_new
    end
end
function close_time_cut(cos)
    global current_sound = current2_sound
    global current_time_range = current2_time_range
    p = Plots.plot(current_time_range,current_sound,xlabel="czas",ylabel = "amplituda")
    Plots.savefig(p, "plot_time.png")
    set_gtk_property!(plot_time_menu,:file,"plot_time.png")
    hide(time_cut_window)
end
function smooth_open(cos)
    set_gtk_property!(smooth_plot,:file,"plot_time.png")
    show(smooth_window)
end
function smooth_sound(cos)
    val = get_gtk_property(smoothness,:value,Int64)
    global current2_sound = NaszModul.smoothing(current_sound,val)
    s = Plots.plot(current_time_range,current2_sound,xlabel="czas",ylabel = "amplituda")
    Plots.savefig(s,"plot_time.png")
    set_gtk_property!(smooth_plot,:file,"plot_time.png")   
end
function smooth_close(cos)
    global current_sound = current2_sound
    set_gtk_property!(plot_time_menu,:file,"plot_time.png") 
    hide(smooth_window)
end
function work_on_time(cos)
    volume = get_gtk_property(adjustment1,:value,Float16)
    speed = get_gtk_property(adjustment2,:value,Float16)
    smooth = get_gtk_property(smoothness,:value,Float64)
    global current2_sound = NaszModul.change_volume(current_sound,volume)
    global current_fs = floor(speed * fs)
    global current_time_range = LinRange(0,length(current_sound)/current_fs,length(current_sound))
    p = Plots.plot(current_time_range,current2_sound,xlabel="czas",ylabel = "amplituda")
    Plots.savefig(p,"plot_time.png")
    set_gtk_property!(main_plot,:file,"plot_time.png")
    set_gtk_property!(plot_time_menu,:file,"plot_time.png")
end

function time_menu_close(cos)
    p = Plots.plot(current_time_range,current_sound,xlabel="czas",ylabel = "amplituda")
    Plots.savefig(p, "plot_time.png")
    set_gtk_property!(main_plot,:file,"plot_time.png")
    global temporary_sound = current2_sound
    global temporary_fs = current_fs
    global temporary_time_range = current2_time_range
    hide(time)
end
#okno częstotliwości
function freq_menu_open(cos)
    global current_sound = temporary_sound
    global current_fs = temporary_fs
    global current2_sound = temporary_sound
    global current_fftw = fft(temporary_sound)
    global temporary_freq = fftfreq(length(current_sound),temporary_fs)
    f = Plots.plot(temporary_freq, abs.(current_fftw),xlims = (0,20000),xlabel="częstotliwości",ylabel = "amplituda")
    Plots.savefig(f,"plot_fftw.png")
    set_gtk_property!(plot_freq_menu,:file,"plot_fftw.png") 
    show(freq_window)
end
function work_on_freq(cos)
    shift = get_gtk_property(freq_shift,:value,Int64)
    global current2_sound =  NaszModul.fft_changer(current_sound,shift)
    current2_freq = fftfreq(length(current_sound), current_fs)
    current2_fftw = fft(current2_sound)
    f = Plots.plot(current2_freq, abs.(current2_fftw),xlims = (0,20000),xlabel="częstotliwości",ylabel = "amplituda")
    Plots.savefig(f,"plot_fftw.png")
    set_gtk_property!(plot_freq_menu,:file,"plot_fftw.png") 

end
function show_filtr_window(cos)
    set_gtk_property!(filtr_plot,:file,"plot_fftw.png")
    show(filtr_window)
end

function high_and_low(cos)
    global current_sound = current2_sound
    highpass_value = get_gtk_property(highpass_scale,:value,Float64)
    lowpass_value = get_gtk_property(lowpass_scale,:value,Float64)
    sound1,fft1, costam, f = NaszModul.lowpass(current_sound,current_fs,lowpass_value)
    global current2_sound = sound1
    sound2 ,fft2,costam2,f2 = NaszModul.highpass(current2_sound,current_fs,highpass_value)
    global current2_sound = sound2
    global current2_fftw = fft(current2_sound)
    current_freq = fftfreq(length(current2_sound),current_fs)
    s = Plots.plot(current_freq, abs.(current2_fftw),xlims = (0,20000),xlabel="częstotliwości",ylabel = "amplituda")
    Plots.savefig(s,"plot_fftw.png")
    set_gtk_property!(filtr_plot,:file,"plot_fftw.png")  
end
function close_filtr_window(cos)
    global current_sound = current2_sound
    global current_fftw = current2_fftw
    set_gtk_property!(plot_freq_menu,:file,"plot_fftw.png")
    hide(filtr_window)
end

function freq_menu_close(cos)
    global temporary_sound = current2_sound
    global temporary_fftw = current2_fftw
    current_time_range = LinRange(0,length(current2_sound)/current_fs,length(current2_sound))
    p = Plots.plot(current_time_range,current2_sound,xlabel="czas",ylabel = "amplituda")
    Plots.savefig(p, "plot_time.png")
    set_gtk_property!(main_plot,:file, "plot_time.png")
    hide(freq_window)
end

function play(cos)
    wavplay(temporary_sound,temporary_fs)
end
function play_current(cos)
    wavplay(current_sound,current_fs)
end
function play_current2(cos)
    wavplay(current2_sound,current_fs)
end

function save(cos)
    file = save_dialog("Save as...",GtkNullContainer() , (GtkFileFilter("*.wav", name="All supported formats"), "*.wav"))
    wavwrite(temporary_sound,temporary_fs, file)

end

b = GtkBuilder(filename="gui\\wersja2.glade")#sciezka do pliku wersja2.glade

win = b["window1"]
btn_time = b["btn_time"]
btn_freq = b["btn_freq"]
btn_file =b["btn_file"]
btn_play = b["play"]
main_plot = b["main_plot"]
time = b["time_window"]
plot_time_menu = b["plot_time_menu"]
volume = b["volume"]
adjustment1=b["adjustment1"]
adjustment2 = b["adjustment2"]
adjustment3 = b["adjustment3"]
adjustment4 = b["adjustment4"]
exit_time = b["exit_time"]
accept_time = b["accept_time"]
play_time = b["play_time"]
freq_window = b["freq_window"]
exit_freq = b["exit_freq"]
accept_freq = b["accept_freq"]
play_freq = b["play_freq"]
plot_freq_menu = b["plot_freq_menu"]
freq_shift = b["adjustment6"]
smoothness = b["smoothness"]
smooth_window = b["smooth_window"]
accept_smooth = b["accept_smooth"]
smooth_button = b["smooth_button"]
highpass_scale = b["highpass_scale"]
lowpass_scale = b["lowpass_scale"]
time_cut_window = b["time_cut_window"]
time_cut_image = b["time_cut_image"]
show_cut_time_btn = b["show_cut_time"]
play_cut_time_btn = b["play_cut_time"]
accept_time_cut = b["accept_time_cut"]
accept_cut_time = b["accept_cut_time"]
filtr_window = b["filtr_window"]
accept_filtr = b["accept_filtr"]
play_filtr = b["play_filtr"]
filtr_back_to_menu = b["filtr_back_to_menu"]
smooth_plot = b["smooth_plot"]
play_smooth =b["play_smooth"]
accept_smooth= b["accept_smooth"]
smooth_back_to_menu = b["smooth_back_to_menu"]
filtr_plot = b["filtr_plot"]
filtr_btn = b["filtr_btn"]
save_to_file = b["save_to_file"]
signal_connect(file_load,btn_file,"clicked")
signal_connect(time_menu_open,btn_time,"clicked")
signal_connect(time_menu_close,exit_time,"clicked")
signal_connect(freq_menu_open,btn_freq,"clicked")
signal_connect(freq_menu_close,exit_freq,"clicked")
signal_connect(play,btn_play,"clicked")
signal_connect(play_current2,play_time,"clicked")
signal_connect(play_current2,play_freq,"clicked")
signal_connect(work_on_time,accept_time,"clicked")
signal_connect(work_on_freq,accept_freq,"clicked")
signal_connect(smooth_open,smooth_button,"clicked")
signal_connect(open_time_cut,show_cut_time_btn,"clicked")
signal_connect(play_current2,play_cut_time_btn,"clicked")
signal_connect(cut_time,accept_time_cut,"clicked")
signal_connect(close_time_cut,accept_cut_time,"clicked")
signal_connect(smooth_sound,accept_smooth,"clicked")
signal_connect(play_current2,play_smooth,"clicked")
signal_connect(smooth_close,smooth_back_to_menu,"clicked")
signal_connect(show_filtr_window,filtr_btn,"clicked")
signal_connect(close_filtr_window,filtr_back_to_menu,"clicked")
signal_connect(play_current2,play_filtr,"clicked")
signal_connect(high_and_low,accept_filtr,"clicked")
signal_connect(save,save_to_file,"clicked")
showall(win)
""