#include "Trail.h"

Trail::Trail(std::istream& fin) : states(), weights()
{
    load(fin);
}

unsigned int Trail::getLength() const
{
    return states.size();
}

void Trail::clear()
{
    states.clear();
    weights.clear();
    totalWeight = 0;
    initialIndex = 0;
}

void Trail::append(const State& state, unsigned int weight)
{
    states.push_back(state);
    weights.push_back(weight);
    totalWeight += weight;
}
/*
void Trail::append(const Trail& otherTrail)
{
    for(unsigned int i=0; i<otherTrail.weights.size(); i++)
        append(otherTrail.states[i], otherTrail.weights[i]);
}
*/
void Trail::prepend(const State& state, unsigned int weight)
{
    states.insert(states.begin(), state);
    weights.insert(weights.begin(), weight);
    totalWeight += weight;
    initialIndex++;
}

void Trail::prepop()
{
    if(initialIndex == 0)
        throw Exception("You probably do not want to prepop the initial two-round trail core...");
    states.erase(states.begin());
    totalWeight -= weights.front();
    weights.erase(weights.begin());
    initialIndex--;
}

void Trail::pop()
{
    states.pop_back();
    totalWeight -= weights.back();
    weights.pop_back();
}

bool Trail::prune(const unsigned int rounds, const unsigned int bestMinimumWeight)
{
    unsigned int sliceWeight = 0;
    unsigned int bestSliceWeight, bestSliceOffset = 0;
    unsigned int i = 0;
    const unsigned int startLength = getLength();

    if(startLength < rounds)
        return false;
    if(startLength == rounds)
        return true;

    for( ; i < rounds; ++i)
        sliceWeight += weights[i];

    bestSliceWeight = sliceWeight;

    for(i = 1 ; i <= startLength - rounds; ++i) {
        sliceWeight -= weights[i - 1];
        sliceWeight += weights[i + rounds - 1];
        if(sliceWeight < bestSliceOffset) {
            bestSliceWeight = sliceWeight;
            bestSliceOffset = i;
        }
    }

    if(bestSliceWeight < bestMinimumWeight) {
        for(i = 0; i < bestSliceOffset; ++i)
            prepop();
        for(i = bestSliceOffset + rounds; i < startLength; ++i)
            pop();
        return true;
    }
    else
        return false;
}

void Trail::set(const Trail& otherTrail)
{
    states = otherTrail.states;
    weights = otherTrail.weights;
    totalWeight = otherTrail.totalWeight;
    initialIndex = otherTrail.initialIndex;
}

std::ostream& operator<<(std::ostream& out, const Trail& trail)
{
    const unsigned int length = trail.getLength();
    if (length == 0) {
        out << "This trail is empty.\n";
        return out;
    }

    out << length << "-round ";
    out << "differential trail core of total weight " << trail.totalWeight << "\n";
    for(int i=COLUMNSIZE-1; i>=0; --i) {
        for(unsigned int j = 0; j<length; ++j) {
            trail.states[j].printRow(out, i);
            out << "    ";
        }
        out << "\n";
    }
    return out;
}

Trail& Trail::operator=(const Trail& t)
{
    if(this != &t) {
        states = t.states;
        weights = t.weights;
        totalWeight = t.totalWeight;
        initialIndex = t.initialIndex;
    }
    return *this;
}

void Trail::save(std::ostream& fout) const
{
    fout << weights.size() << ' ';
    fout << totalWeight << ' ';
    fout << initialIndex << ' ';
    for(unsigned int i=0; i<weights.size(); i++)
        fout << weights[i] << ' ';
    fout << states.size() << ' ';
    for(unsigned int i=0; i<states.size(); i++)
        for(unsigned int j=0;j<COLUMNSIZE; j++)
            fout << states[i][j] << ' ';
    fout << '\n';
}

void Trail::load(std::istream& fin)
{
    unsigned int size;
    fin >> size;
    if ((size == 0) || (fin.eof()))
        throw TrailException(std::string("Could not read file"));
    fin >> totalWeight;

    fin >> initialIndex;

    weights.resize(size);
    for(unsigned int i=0; i<size; i++)
        fin >> weights[i];

    states.resize(size);
    for(unsigned int i=0; i<size; i++) {
        State s;
        for(unsigned int j=0; j<COLUMNSIZE; j++)
            fin >> s[j];
        states[i] = s;
    }
    fin.ignore(10,'\n');
}

bool Trail::loadBruteforce(std::istream& fin) {
    char buf[16];
    State a, b;

    fin.read(buf, 16);

    if(!fin.good()) {
        return false;
    }
    for(unsigned int i = 0; i < 8; ++i) {
        a.setColumn(buf[i] & 0xf, 2*i);
        a.setColumn(buf[i] >> 4, 2*i+1);
    }
    for(unsigned int i = 0; i < 8; ++i) {
        b.setColumn(buf[8+i] & 0xf, 2*i);
        b.setColumn(buf[8+i] >> 4, 2*i+1);
    }

    clear();
    applyDispersion(b);
    append(a, a.hammingWeight());
    append(b, b.hammingWeight());
    return true;
}
