#ifndef MUSFUN_A1812399_2ECF_4992_A7ED_983C8B03FD54
#define MUSFUN_A1812399_2ECF_4992_A7ED_983C8B03FD54

#include <vector>
#include <string>
#include <muParser.h>
#include <iostream>
#include <limits>
#include <concepts>

/**
 * @brief Class that implements muparser-interfaced functions R^ --> R
 */

template <std::floating_point T>
class mu_Sfun
{

private:
    std::array<T, 2> m_vars;
    mu::Parser m_parser;

public:
    // Constructor
    mu_Sfun(const std::string &s)
    {
        try
        {
            m_parser.DefineVar("x1", &m_vars[0]);
            m_parser.DefineVar("x2", &m_vars[1]);
            m_parser.SetExpr(s);
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cerr << e.GetMsg() << std::endl;
        }
    }

    // Copy constructor
    mu_Sfun(const mu_Sfun &other)
        : m_vars(other.m_vars), m_parser(other.m_parser)
    {
        m_parser.DefineVar("x1", &m_vars[0]);
        m_parser.DefineVar("x2", &m_vars[1]);
    }

    // Move constructor
    mu_Sfun(mu_Sfun &&other) noexcept
        : m_vars(std::move(other.m_vars)), m_parser(std::move(other.m_parser))
    {
        m_parser.DefineVar("x1", &m_vars[0]);
        m_parser.DefineVar("x2", &m_vars[1]);
    }

    // Copy assignment operator
    mu_Sfun &operator=(const mu_Sfun &other)
    {
        if (this != &other)
        {
            m_vars = other.m_vars;
            m_parser = other.m_parser;
            m_parser.DefineVar("x1", &m_vars[0]);
            m_parser.DefineVar("x2", &m_vars[1]);
        }
        return *this;
    }

    // Move assignment operator
    mu_Sfun &operator=(mu_Sfun &&other) noexcept
    {
        if (this != &other)
        {
            m_vars = std::move(other.m_vars);
            m_parser = std::move(other.m_parser);
            m_parser.DefineVar("x1", &m_vars[0]);
            m_parser.DefineVar("x2", &m_vars[1]);
        }
        return *this;
    }

    // Destructor
    ~mu_Sfun() = default;

    T operator()(const std::array<T, 2> &x)
    {
        m_vars = x;
        T y = m_parser.Eval();
        return y;
    }
};

#endif /* MUXFUNCT_A1812399_2ECF_4992_A7ED_983C8B03FD54 */
