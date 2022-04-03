template<typename T> static void read(T& var) { std::cin >> var; }
template<typename T, typename... Ts> static void read(T& var, Ts&... vars) { read(var); read(vars...); }
