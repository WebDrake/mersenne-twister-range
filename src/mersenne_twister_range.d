/++
The Mersenne Twister generator.

Copyright: Copyright Andrei Alexandrescu 2008 - 2009, Ilya Yaroshenko 2016-.
License:   $(HTTP www.boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: $(HTTP erdani.org, Andrei Alexandrescu) Ilya Yaroshenko (rework)
+/
module mersenne_twister_range;

// Segments of the code in this file Copyright (c) 1997 by Rick Booth
// From "Inner Loops" by Rick Booth, Addison-Wesley

// Work derived from:

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

import std.range.primitives;
import std.traits;

/++
The $(LUCKY Mersenne Twister) generator.
+/
struct MersenneTwisterEngine(UIntType, size_t w, size_t n, size_t m, size_t r,
                             UIntType a, size_t u, UIntType d, size_t s,
                             UIntType b, size_t t,
                             UIntType c, size_t l, UIntType f)
    if (isUnsigned!UIntType)
{
    static assert(0 < w && w <= UIntType.sizeof * 8);
    static assert(1 <= m && m <= n);
    static assert(0 <= r && 0 <= u && 0 <= s && 0 <= t && 0 <= l);
    static assert(r <= w && u <= w && s <= w && t <= w && l <= w);
    static assert(0 <= a && 0 <= b && 0 <= c);
    static assert(n <= sizediff_t.max);

    ///Mark this as a Rng
    enum bool isUniformRandom = true;

/**
Parameters for the generator.
*/
    enum size_t   wordSize   = w;
    enum size_t   stateSize  = n; /// ditto
    enum size_t   shiftSize  = m; /// ditto
    enum size_t   maskBits   = r; /// ditto
    enum UIntType xorMask    = a; /// ditto
    enum size_t   temperingU = u; /// ditto
    enum UIntType temperingD = d; /// ditto
    enum size_t   temperingS = s; /// ditto
    enum UIntType temperingB = b; /// ditto
    enum size_t   temperingT = t; /// ditto
    enum UIntType temperingC = c; /// ditto
    enum size_t   temperingL = l; /// ditto
    enum UIntType initializationMultiplier = f; /// ditto

    /// Smallest generated value (0).
    enum UIntType min = 0;
    /// Largest generated value.
    enum UIntType max = UIntType.max >> (UIntType.sizeof * 8u - w);
    // note, `max` also serves as a bitmask for the lowest `w` bits
    static assert(a <= max && b <= max && c <= max && f <= max);

    /// The default seed value.
    enum UIntType defaultSeed = 5489u;

    /// Bitmasks used in the 'twist' part of the algorithm
    private enum UIntType lowerMask = (cast(UIntType) 1u << r) - 1;
    private enum UIntType upperMask = (~lowerMask) & this.max; // ditto

    /// Collection of all state variables
    /// used by the generator
    private struct State
    {
        /// State array of the generator
        UIntType[n] data;

        /// Cached copy of most recently updated
        /// element of `data` state array, ready
        /// to be tempered to generate next
        /// `front` value
        UIntType z;

        /// Most recently generated random variate
        UIntType front;

        /// Index of the entry in the `data`
        /// state array that will be twisted
        /// in the next `popFront()` call
        size_t index;
    }

    /// State variables used by the generator;
    /// initialized to values equivalent to
    /// explicitly seeding the generator with
    /// `defaultSeed`
    private State state = defaultState();
    // NOTE: the above is a workaround to ensure
    // backwards compatibility with the original
    // implementation, which permitted implicit
    // construction.  With `@disable this();`
    // it would not be necessary.

/**
   Constructs a MersenneTwisterEngine object.
*/
    this(UIntType value) @safe pure nothrow @nogc
    {
        seed(value);
    }

    /**
       Generates the default initial state for a Mersenne
       Twister; equivalent to the internal state obtained
       by calling `seed(defaultSeed)`
    */
    private static State defaultState() @safe pure nothrow @nogc
    {
        if (!__ctfe) assert(false);
        State mtState;
        seedImpl(defaultSeed, mtState);
        return mtState;
    }

/**
   Seeds a MersenneTwisterEngine object.
   Note:
   This seed function gives 2^w starting points (the lowest w bits of
   the value provided will be used). To allow the RNG to be started
   in any one of its internal states use the seed overload taking an
   InputRange.
*/
    void seed()(UIntType value = defaultSeed) @safe pure nothrow @nogc
    {
        this.seedImpl(value, this.state);
    }

    /**
       Implementation of the seeding mechanism, which
       can be used with an arbitrary `State` instance
    */
    private static void seedImpl(UIntType value, ref State mtState)
    {
        mtState.data[$-1] = value;
        static if (this.max != UIntType.max)
        {
            mtState.data[$-1] &= this.max;
        }

        foreach_reverse (size_t i, ref e; mtState.data[0 .. $-1])
        {
            e = f * (mtState.data[i + 1] ^ (mtState.data[i + 1] >> (w - 2))) + cast(UIntType)(n - (i + 1));
            static if (this.max != UIntType.max)
            {
                e &= this.max;
            }
        }

        mtState.index = n-1;

        // double popFront() to guarantee both `mtState.z`
        // and `mtState.front` are derived from the newly
        // set values in `mtState.data`
        MersenneTwisterEngine.popFrontImpl(mtState);
        MersenneTwisterEngine.popFrontImpl(mtState);
    }

/**
   Seeds a MersenneTwisterEngine object using an InputRange.

   Throws:
   $(D Exception) if the InputRange didn't provide enough elements to seed the generator.
   The number of elements required is the 'n' template parameter of the MersenneTwisterEngine struct.
 */
    void seed(T)(T range) if (isInputRange!T && is(Unqual!(ElementType!T) == UIntType))
    {
        this.seedImpl(range, this.state);
    }

    /**
       Implementation of the range-based seeding mechanism,
       which can be used with an arbitrary `State` instance
    */
    private static void seedImpl(T)(T range, ref State mtState)
        if (isInputRange!T && is(Unqual!(ElementType!T) == UIntType))
    {
        size_t j;
        for (j = 0; j < n && !range.empty; ++j, range.popFront())
        {
            sizediff_t idx = n - j - 1;
            mtState.data[idx] = range.front;
        }

        mtState.index = n - 1;

        if (range.empty && j < n)
        {
            import core.internal.string : UnsignedStringBuf, unsignedToTempString;

            UnsignedStringBuf buf = void;
            string s = "MersenneTwisterEngine.seed: Input range didn't provide enough elements: Need ";
            s ~= unsignedToTempString(n, buf, 10) ~ " elements.";
            throw new Exception(s);
        }

        // double popFront() to guarantee both `mtState.z`
        // and `mtState.front` are derived from the newly
        // set values in `mtState.data`
        MersenneTwisterEngine.popFrontImpl(mtState);
        MersenneTwisterEngine.popFrontImpl(mtState);
    }

/**
   Advances the generator.
*/
    void popFront() @safe pure nothrow @nogc
    {
        this.popFrontImpl(this.state);
    }

    /// Internal implementation of `popFront()`, which
    /// can be used with an arbitrary `State` instance
    private static void popFrontImpl(ref State mtState)
    {
        // This function blends two nominally independent
        // processes: (i) calculation of the next random
        // variate `mtState.front` from the cached previous
        // `data` entry `z`, and (ii) updating the value
        // of `data[index]` and `mtState.z` and advancing
        // the `index` value to the next in sequence.
        //
        // By interweaving the steps involved in these
        // procedures, rather than performing each of
        // them separately in sequence, the variables
        // are kept 'hot' in CPU registers, allowing
        // for significantly faster performance.
        sizediff_t index = mtState.index;
        sizediff_t next = index - 1;
        if (next < 0)
            next = n - 1;
        auto z = mtState.z;
        sizediff_t conj = index - m;
        if (conj < 0)
            conj = index - m + n;

        static if (d == UIntType.max)
        {
            z ^= (z >> u);
        }
        else
        {
            z ^= (z >> u) & d;
        }

        auto q = mtState.data[index] & upperMask;
        auto p = mtState.data[next] & lowerMask;
        z ^= (z << s) & b;
        auto y = q | p;
        auto x = y >> 1;
        z ^= (z << t) & c;
        if (y & 1)
            x ^= a;
        auto e = mtState.data[conj] ^ x;
        z ^= (z >> l);
        mtState.z = mtState.data[index] = e;
        mtState.index = next;

        // technically we should take the lowest `w`
        // bits here, but if the tempering bitmasks
        // `b` and `c` are set correctly, this should
        // be unnecessary
        mtState.front = z;
    }

/**
   Returns the current random value.
 */
    @property UIntType front() @safe const pure nothrow @nogc
    {
        return this.state.front;
    }

///
    @property typeof(this) save() @safe pure nothrow @nogc
    {
        return this;
    }

/**
Always $(D false).
 */
    enum bool empty = false;
}

/**
A $(D MersenneTwisterEngine) instantiated with the parameters of the
original engine $(HTTP math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html,
MT19937), generating uniformly-distributed 32-bit numbers with a
period of 2 to the power of 19937. Recommended for random number
generation unless memory is severely restricted, in which case a $(D
LinearCongruentialEngine) would be the generator of choice.
 */
alias Mt19937 = MersenneTwisterEngine!(uint, 32, 624, 397, 31,
                                       0x9908b0df, 11, 0xffffffff, 7,
                                       0x9d2c5680, 15,
                                       0xefc60000, 18, 1812433253);

///
@safe unittest
{
    import std.random : unpredictableSeed;
    // seed with a constant
    Mt19937 gen;
    auto n = gen.front; // same for each run
    // Seed with an unpredictable value
    gen.seed(unpredictableSeed);
    n = gen.front; // different across runs
}

@safe nothrow unittest
{
    import std.algorithm;
    import std.random : isUniformRNG, isSeedable, unpredictableSeed;
    import std.range;
    static assert(isUniformRNG!Mt19937);
    static assert(isUniformRNG!(Mt19937, uint));
    static assert(isSeedable!Mt19937);
    static assert(isSeedable!(Mt19937, uint));
    static assert(isSeedable!(Mt19937, typeof(map!((a) => unpredictableSeed)(repeat(0)))));
    Mt19937 gen;
    assert(gen.front == 3499211612);
    popFrontN(gen, 9999);
    assert(gen.front == 4123659995);
    try { gen.seed(iota(624u)); } catch (Exception) { assert(false); }
    assert(gen.front == 3708921088u);
    popFrontN(gen, 9999);
    assert(gen.front == 165737292u);
}

/**
A $(D MersenneTwisterEngine) instantiated with the parameters of the
original engine $(HTTP en.wikipedia.org/wiki/Mersenne_Twister,
MT19937-64), generating uniformly-distributed 64-bit numbers with a
period of 2 to the power of 19937.
*/
alias Mt19937_64 = MersenneTwisterEngine!(ulong, 64, 312, 156, 31,
                                          0xb5026f5aa96619e9, 29, 0x5555555555555555, 17,
                                          0x71d67fffeda60000, 37,
                                          0xfff7eee000000000, 43, 6364136223846793005);

///
@safe unittest
{
    import std.random : unpredictableSeed;
    // Seed with a constant
    auto gen = Mt19937_64(12345);
    auto n = gen.front; // same for each run
    // Seed with an unpredictable value
    gen.seed(unpredictableSeed);
    n = gen.front; // different across runs
}

@safe nothrow unittest
{
    import std.algorithm;
    import std.random : isUniformRNG, isSeedable, unpredictableSeed;
    import std.range;
    static assert(isUniformRNG!Mt19937_64);
    static assert(isUniformRNG!(Mt19937_64, ulong));
    static assert(isSeedable!Mt19937_64);
    static assert(isSeedable!(Mt19937_64, ulong));
    // FIXME: this test demonstrates viably that Mt19937_64
    // is seedable with an infinite range of `ulong` values
    // but it's a poor example of how to actually seed the
    // generator, since it can't cover the full range of
    // possible seed values.  Ideally we need a 64-bit
    // unpredictable seed to complement the 32-bit one!
    static assert(isSeedable!(Mt19937_64, typeof(map!((a) => (cast(ulong)unpredictableSeed))(repeat(0)))));
    Mt19937_64 gen;
    assert(gen.front == 14514284786278117030uL);
    popFrontN(gen, 9999);
    assert(gen.front == 9981545732273789042uL);
    try { gen.seed(iota(312uL)); } catch (Exception) { assert(false); }
    assert(gen.front == 14660652410669508483uL);
    popFrontN(gen, 9999);
    assert(gen.front == 15956361063660440239uL);
}
