"""Funkcja zwacająca macierz DFT"""
function dft(signal)
    N = length(signal)
    zeta_powers = OffsetArray([exp( -2π *  im * n / N) for n in 0:(N-1) ], 0:(N-1))
    [  sum( signal[n + 1] * zeta_powers[(n * f) % N] for n in 0:(N-1)   ) for f in 0:(N-1) ]
end

"""Funkcja zwracająca częstotliwość próbkowania DFT"""
function dftfreq(n, fs)
    if n >= 1
        f = zeros(n)
        if iseven(n)
            positive = collect(0:((n/2)-1))
            negative = collect((-n/2):1:(-1))
        else isodd(N)
            positive = collect(0:((n-1)/2))
            negative = collect((-(n-1)/2):1:(-1))
        end
        len = length(positive)
        f[1:len] = positive
        f[(len + 1):end] = negative
        return f*fs/n
    else
        throw(ArgumentError("n powinno być większe równe 1"))
    end
end

"""Odwrotna funkcja DFT"""
function idft(signal)
    N = length(signal)
    zeta_powers = OffsetArray([exp( 2π *  im * n / N) for n in 0:(N-1) ], 0:(N-1))
    [  (1/N)*sum( signal[n + 1] * zeta_powers[(n * f) % N] for n in 0:(N-1)   ) for f in 0:(N-1) ]
end

"""Funkcja wycinająca mniejsze wielkości"""
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

"""Funkcja wycinająca większe wielkości"""
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

"""Funkcja wycinająca określoną częstotliwość"""
function remove_frequency(sound, freq_to_remove, dt = 0.001, start = 0, stop = 220.5)
    t = LinRange(start, stop, Int((stop - start)/dt))
    n = length(t) 
    fhat = real.(fft(sound[1][:, 1])) 
    PSD = fhat .* conj(fhat)/n 
    freq = 1/(dt * n) .* [i for i in 1:n] 
    
    indices = [i*dt*n != freq_to_remove for i in freq] 
    PSD_clean = indices .* PSD 
    fhat_clean = indices .* fhat 
    real(ifft(fhat_clean)), freq, PSD_clean, fhat_clean
end

"""Funkcja do przyspieszania(?)"""
function change_speed(sound, fs, speed, file_name)
    wavwrite(sound, file_name, Fs = fs * speed)
    sound_speeded, fs_speeded = wavread(file_name)
    return sound_speeded, fs_speeded
end

"""Funkcja do zmiany głośniści(?)"""
function change_volume(sound, volume)
    sound_fft = fft(sound)
    new_sound_fft = volume .* sound_fft
    new_sound = ifft(new_sound_fft)
    return new_sound
end
