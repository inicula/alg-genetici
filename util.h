template<typename T>
static void read(T& var)
{
        finput >> var;
}

template<typename T, typename... Ts>
static void read(T& var, Ts&... vars)
{
        read(var);
        read(vars...);
}

template<typename Cont, typename Key>
bool contains(const Cont& c, const Key& k)
{
        return c.find(k) != c.cend();
}
