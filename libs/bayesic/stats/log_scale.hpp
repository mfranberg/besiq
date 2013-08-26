#ifndef __LOG_SCALE_H__
#define __LOG_SCALE_H__

#include <cmath>
#include <cfloat>
#include <stdexcept>

template <class T>
class log_scale
{
public:
    log_scale(T value)
    {
        if( value != 0.0 )
        {
            m_log_value = log( value );
            m_is_zero = false;
        }
        else
        {
            m_log_value = -DBL_MAX;
            m_is_zero = true;
        }
    }

    log_scale(const log_scale<T> &other)
    {
        m_log_value = other.m_log_value;
        m_is_zero = other.m_is_zero;
    }

    template <class K>
    static log_scale<K> from_log(K log_value)
    {
        log_scale<K> from_log( 1.0 );
        from_log.m_log_value = log_value;

        return from_log;
    }

    const log_scale<T> &operator=(const log_scale<T> &other)
    {
        if( this != &other )
        {
            m_log_value = other.m_log_value;
            m_is_zero = other.m_is_zero;
        }
        return *this;
    }

    const log_scale<T> &operator+=(const log_scale<T> &other)
    {
        if( !m_is_zero && !other.m_is_zero )
        {
            T diff = other.m_log_value - m_log_value;
            m_log_value += log( 1.0 + exp( diff ) );
        }
        else if( m_is_zero )
        {
            m_log_value = other.m_log_value;
            m_is_zero = other.m_is_zero;
        }

        return *this;
    }

    const log_scale<T> &operator-=(const log_scale<T> &other)
    {
        if( !m_is_zero && !other.m_is_zero )
        {
            T diff = m_log_value - other.m_log_value;
            if( diff < 0.0 )
            {
                throw std::domain_error( "log_scale: Cannot represent negative values on log scale." );
            }
            m_log_value = other.m_log_value + log( exp( diff ) - 1.0 );
        }
        else if( m_is_zero && !other.m_is_zero )
        {
            throw std::domain_error( "log_scale: Cannot represent negative values on log scale." );
        }

        return *this;
    }

    const log_scale<T> &operator*=(const log_scale<T> &other)
    {
        if( !m_is_zero && !other.m_is_zero )
        {
            m_log_value += other.m_log_value;
        }
        else
        {
            m_log_value = -DBL_MAX;
            m_is_zero = true;
        }

        return *this;
    }

    const log_scale<T> &operator/=(const log_scale<T> &other)
    {
        /* Note: Order matters */
        if( !m_is_zero && !other.m_is_zero )
        {
            m_log_value -= other.m_log_value;
        }
        else if( other.m_is_zero )
        {
            throw std::domain_error( "log_scale: Division with zero." );
        }
        else
        {
            m_log_value = -DBL_MAX;
            m_is_zero = true;
        }

        return *this;
    }

    T value()
    {
        if( !m_is_zero )
        {
            return exp( m_log_value );
        }
        else
        {
            return 0.0;
        }
    }

    T log_value()
    {
        if( !m_is_zero )
        {
            return m_log_value;
        }
        else
        {
            return -DBL_MAX;
        }
    }

private:
    /**
     * The underlying logged value.
     */
    T m_log_value;

    /**
     * Determines whether this should be zero.
     */
    bool m_is_zero;
};

template <class T>
log_scale<T> operator+(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) += rhs;
}

template <class T>
log_scale<T> operator+(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) += rhs;
}

template <class T>
log_scale<T> operator+(const log_scale<T> &lhs, const T &rhs)
{
    return log_scale<T>( lhs ) += log_scale<T>( rhs );
}

template <class T>
log_scale<T> operator-(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) -= rhs;
}

template <class T>
log_scale<T> operator-(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) -= rhs;
}

template <class T>
log_scale<T> operator-(const log_scale<T> &lhs, const T &rhs)
{
    return log_scale<T>( lhs ) -= log_scale<T>( rhs );
}


template <class T>
log_scale<T> operator/(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) /= rhs;
}

template <class T>
log_scale<T> operator/(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) /= rhs;
}

template <class T>
log_scale<T> operator/(const log_scale<T> &lhs, const T &rhs)
{
    return log_scale<T>( lhs ) /= log_scale<T>( rhs );
}

template <class T>
log_scale<T> operator*(const log_scale<T> &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) *= rhs;
}

template <class T>
log_scale<T> operator*(const T &lhs, const log_scale<T> &rhs)
{
    return log_scale<T>( lhs ) *= rhs;
}

template <class T>
log_scale<T> operator*(const log_scale<T> &lhs, const T &rhs)
{
    return log_scale<T>( lhs ) *= log_scale<T>( rhs );
}

typedef log_scale<double> log_double;
typedef log_scale<float> log_float;
typedef log_scale<int> log_int;

#endif /* End of __LOG_SCALE_H__ */
