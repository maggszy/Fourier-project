"""
Informacje:
aby użyć modułu piszemy include("naszmodul.jl") (tak jak te wszystkie usingi)
A później używamy tych funkcji robiąc np. NaszModul.dftfreq(coś tam coś tam)
Nazwa NaszModul jest oczywiście dość robocza...
"""

"""Moduł z funkcjami do przetwarzania dźwięku"""
module NaszModul

using OffsetArrays 
using WAV
using FFTW

"""
    dft(signal)

Funkcja zwacająca tablice DFT(rozkład amplitudy od  częstotliwości), gdzie 'signal'
to tablica zawierająca wartości amplitudy w zależności od czasu. 
"""
function dft(signal)
    N = length(signal)
    zeta_powers = OffsetArray([exp( -2π *  im * n / N) for n in 0:(N-1) ], 0:(N-1))
    [  sum( signal[n + 1] * zeta_powers[(n * f) % N] for n in 0:(N-1)   ) for f in 0:(N-1) ]
end

"""
    dftfreq(n,fs)

Funkcja zwracająca częstotliwość próbkowania DFT, gdzie 'n' to ilość próbek,
a 'fs' częśtotliwość próbkowania sygnału wejściowego.
"""
function dftfreq(n::Int64, fs::Float64)
    if n >= 1
        if iseven(n)
            positive = 0:((n/2)-1)
            negative = (-n/2):(-1)
        else
            positive = 0:((n-1)/2)
            negative = (-(n-1)/2):(-1)
        end
        f = [positive; negative]
        return f*fs/n
    else
        throw(ArgumentError("n powinno być większe równe 1"))
    end
end

"""
    idft(signal)

Funkcja zwracająca odrtotne DFT(zmienia rozkład zależny od częstotliwości
na taki zależny od czasu), gdzie 'signal' to tablica zawierająca zależności
amplitudy od częstotliwości.
"""
function idft(signal)
    N = length(signal)
    zeta_powers = OffsetArray([exp( 2π *  im * n / N) for n in 0:(N-1) ], 0:(N-1))
    [  (1/N)*sum( signal[n + 1] * zeta_powers[(n * f) % N] for n in 0:(N-1)   ) for f in 0:(N-1) ]
end

function fft2(x)
    N = length(x)
    if log2(N) % 1 > 0
        throw(ArgumentError("must be a power of 2"))
    end
    N_min = min(N, 32)
    n = 0:(N_min-1)
    M = cis.(-2π * n * n'/N_min)
    X = M * reshape(x, (:, N_min))'
    while size(X)[1] < N
        X_even = X[:, 1:Int(size(X)[2] / 2)]
        X_odd = X[:, (Int(size(X)[2] / 2)+1):end]
        terms = cis.(-1 * π * collect(0:(size(X)[1]-1))/size(X)[1])
        X = [X_even .+ terms .* X_odd; X_even .- terms .* X_odd]
    end
    return vec(X)
end

function ifft2(x)
    N = length(x)
    if log2(N) % 1 > 0
        throw(ArgumentError("must be a power of 2"))
    end
    N_min = min(N, 32)
    
    n = 0:(N_min-1)
    M = cis.(2π * n * n'/ N_min)
    X = M * reshape(x, (:, N_min))'
    while size(X)[1] < N
        X_even = X[:, 1:Int(size(X)[2] / 2)]
        X_odd = X[:, (Int(size(X)[2] / 2)+1):end]
        terms = cis.(π * collect(0:(size(X)[1]-1))/size(X)[1])
        X = [X_even .+ terms .* X_odd; X_even .- terms .* X_odd]
    end
    return (1/N).*[X[1]; reverse(X[2:end])]
end

"""
    pow2matrix(x)

Funkcja rozszerzająca macierz o zera
tak, by długość była potęgą liczby 2"""
function pow2matrix(x)
    N = length(x)
    n = nextpow(2, N)
    M = zeros(n)
    M[1:N] = x
    return M, N
end

"""
    fft_w(x)

Funkcja używająca algorytmu FFT, by zwrócić macierz rozkładu amplitudy
od częstotliwości dla macierzy "x" zawierającej wartości amplitudy w zależności od czasu. 
"""
function fft_w(x)
    pow2_signal = pow2matrix(x, Float64)
    x = pow2_signal[1]
    N = length(x)
    X = fft2(x)
    return vec(X), pow2_signal[2]
end

"""
    freqincrease(signal::Array{Float64,1}, shift::Int64)

Funkcja zwiększa częstotliwość tablicy zawierającej wartości amplitudy w zależności
od czasu (signal) o podaną liczbę cykli (shift).
"""
function freqincrease(signal::Array{Float64,1}, shift::Int64)
    fft_signal = fft(signal)
    amplitude = abs.(fft_signal)
    angle_signal = angle.(fft_signal)
    N = length(fft_signal)
    amp1 = amplitude[2:Int(N/2)+1]                
    angle1 = angle_signal[2:Int(N/2)+1] 
    amp2 = amplitude[Int(N/2)+2:end]
    angle2 = angle_signal[Int(N/2)+2:end] 
    amp1shift   = [zeros(shift); amp1[1:end-shift]]
    angle1shift = [zeros(shift); angle1[1:end-shift]]
    amp2shift   = [amp2[shift+1:end]; zeros(shift)]
    angle2shift = [angle2[shift+1:end]; zeros(shift)]
    amplitude_new = [amplitude[1]; amp1shift; amp2shift]
    angle_new  = [angle_signal[1]; angle1shift; angle2shift]
    x = amplitude_new .* cos.(angle_new)
    y = amplitude_new .* sin.(angle_new)
    fft_new = x + im*y
    return ifft(fft_new)
end

"""
    freqreduce(signal::Array{Float64,1}, shift::Int64)

Funkcja zwiększa częstotliwość tablicy zawierającej wartości amplitudy w zależności
od czasu (signal) o podaną liczbę cykli (shift).
"""
function freqreduce(signal::Array{Float64,1}, shift::Int64)
    fft_signal = fft(signal)
    amplitude = abs.(fft_signal)
    angle_signal = angle.(fft_signal)
    N = length(fft_signal)
    amp1 = amplitude[2:Int(N/2)+1]                
    angle1 = angle_signal[2:Int(N/2)+1] 
    amp2 = amplitude[Int(N/2)+2:end]
    angle2 = angle_signal[Int(N/2)+2:end]
    amp1shift   = [amp1[shift+1:end]; zeros(shift)]
    angle1shift = [angle1[shift+1:end]; zeros(shift)]
    amp2shift   = [zeros(shift); amp2[1:end-shift]]
    angle2shift = [zeros(shift); angle2[1:end-shift]]
    amplitude_new = [amplitude[1]; amp1shift; amp2shift]
    angle_new  = [angle_signal[1]; angle1shift; angle2shift]
    x = amplitude_new .* cos.(angle_new)
    y = amplitude_new .* sin.(angle_new)
    fft_new = x + im*y
    return ifft(fft_new)
end

function fft_changer(signal::Array{Float64,1}, shift::Int64)
    if shift >= 0
        freqincrease(signal, shift)
    else
        freqreduce(signal, abs(shift))
    end
end
"""
    denoising_high(sound,fs,N,cut_point)

Funkcja wycinająca mniejsze wielkości(tzn. cichsze, na podstawie PSD)
"""
function denoising_high(sound, fs, N,cut_point)
    tstep=1/fs
    time=LinRange(0, (N-1)/fs, N)
    freq_step = fs/N
    freq = LinRange(0, (N-1)*freq_step, N)
    
    fhat=fft(sound)
    PSD = real.(fhat .* conj(fhat)/N) #power spectral density
    
    indices=[i > cut_point for i in PSD]
    PSD_clean= indices .* PSD
    fhat_clean = indices .* fhat
    
    return time, PSD_clean, fhat_clean, freq
end

"""
    denoising_low(sound,fs,N,cut_point)

Funkcja wycinająca większe wielkości(tzn. głośniejsze, na podstawie PSD)
"""
function denoising_low(sound, fs, N,cut_point)
    tstep=1/fs
    time=LinRange(0, (N-1)/fs, N) 
    freq_step = fs/N
    freq = LinRange(0, (N-1)/fs, N)
    
    fhat=fft(sound)
    PSD = real.(fhat .* conj(fhat)/N)
    
    indices=[i < cut_point for i in PSD]
    PSD_clean= indices .* PSD
    fhat_clean = indices .* fhat

    return time, PSD_clean, fhat_clean, freq
end

"""
    remove_frequency(sound,fs,N,freq_start,freq_stop)

Funkcja wycinająca określony zakres częstotliwości
"""
function remove_frequency(sound, fs, N, freq_start, freq_stop)
    tstep = 1/fs # sample time interval
    t = LinRange(0, (N-1)*tstep, N) # time steps
    fstep = fs/N # freq interval
    f = LinRange(0, (N-1)*fstep, N) # freq steps
    
    fhat = fft(sound)
    fhat[floor(Int, freq_start * N/f[end]):ceil(Int, freq_stop * N/f[end])] .= 0
    fhat[floor(Int, (f[end] - freq_stop) * N/f[end]):ceil(Int, (f[end] - freq_start) * N/f[end])] .= 0

    fhat_mag = abs.(fhat)/N
    return fhat, fhat_mag, f
end

"""
    change_speed(sound,fs,speed,file_name)

Funkcja do przyspieszania/zwalniania tempa odtwarzania, gdzie 'sound'
to tablica zawierająca rozkład amplitudy w zależności od czasu, fs to
częstotliwość próbkowania sygnału wejściowego, 'speed' to współczynnik 
prędkości dźwięu w nowym pliku 
"""
function change_speed(sound, fs, speed, file_name)
    wavwrite(sound, file_name, Fs = fs * speed)
    sound_speeded, fs_speeded = wavread(file_name)
    return sound_speeded, fs_speeded
end

"""
    change_volume(sound,volume)

Funkcja do zmiany głośniści (ściszanie i zgłaśnianie), gdzie 'sound'
to tablica zawierająca rozkład amplitudy w zależności od czasu, a volume
to współczynnnik głośności dźwięku w nowym pliku. 
"""
function change_volume(sound, volume)
    new_sound=volume*sound
    return new_sound
end

"""
    cutting_time(sound,fs,start,stop)

Funkcja do przycinania czasu nagrania, gdzie 'sound'to tablica zawierająca
rozkład amplitudy w zależności od czasu, fs to częstotliwość próbkowania 
sygnału wejściowego, a start i stop to miejsca w czasie podane w sekundach
"""
function cutting_time(sound,fs,start=0,stop=((length(sound[:,1])-1)/fs))
    N=length(sound[:,1]) 
    time=LinRange(0,(N-1)/fs,N)
    
    if start==0
        new_sound=sound[:,1][1:floor(Int,stop/time[end]*N)]
    else 
        new_sound=sound[:,1][floor(Int,start/time[end]*N):floor(Int,stop/time[end]*N)]
    end
        
    N_new=length(new_sound)
    time_new=LinRange(0,(N_new-1)/fs,N_new)
    
    return time_new, new_sound
end

"""
    loess(index,X,Y,m)

Oszacowanie wartości w zadanym (wartością 'index') punkcie z tablicy 'X',
za pomocą regresji lokalnej. Y to zbiór wartości odpowiadających tablicy'X', a 2m+1 to
szerokośc otoczenia
"""
function loess(index, X,Y,m)
    if index <= m
        xs = X[1:(index+m)]
        ys = Y[1:(index+m)]
    elseif index + m <= length(X)
        xs = X[(index-m):(index+m)]
        ys = Y[(index-m):(index+m)]
    else
        xs = X[(index-m):end]
        ys = Y[(index-m):end]    
    end
    av_x = sum(xs)/length(xs) 
    av_y = sum(ys)/length(ys)
    
    a = sum((xs .- av_x) .* (ys .- av_y) )/ sum( (xs .- av_x).^2)
    b = av_y - a * av_x
    return a*X[index]+b
end

"""
    smoothing(Y, m)

Funkcja, która za pomocą loess znajduje przybliżenie funkcji zadanej przez tablicę 'Y'.
Uzywamy do odszumiania nagran.
"""
function smoothing(Y,m= 10)
    X = 1:length(Y)
    y = [loess(index, X,Y,m) for index in 1:length(X)]
    return y
end

end