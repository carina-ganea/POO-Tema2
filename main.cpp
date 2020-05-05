#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <exception>
#include <cmath>
#include <cfloat>
#include <iomanip>
using namespace std;

#define MAX 10000
#define MAX_CIUR 200
#define MAX_ITER 500

bool ciur[MAX_CIUR];
void ciur_eratostene(){
    for( int i = 2; i < MAX_CIUR; i++){
        int pas = i * 2;
        while( pas < MAX_CIUR){
            ciur[pas] = true;
            pas +=i;
        }
    }
}
class ExceptieEisenstein: public exception
{
public:
    virtual const char* what() const throw()
    {
        return "Criteriul lui Eisenstein nu poate fi aplicat, prea multe zecimale\n";
    }

} excpE;


class Monom{
    int grad_;
    float coef_;
public:
    Monom(): grad_(0), coef_(1){}
    explicit Monom(const int& Grad, const float& Coef = 1): grad_(Grad), coef_(Coef){}
    Monom( const Monom& m)  = default;
    ~Monom()= default;
    Monom& operator= ( const Monom& m){
        grad_ = m.grad();
        coef_ = m.coef();
        return *this;
    }
    int grad() const{
        return grad_;
    }
    float coef() const{
        return coef_;
    }
    void setCoef( float x){
        coef_ = x;
    }
    void setGrad( int x){
        grad_ = x;
    }
};

istream& operator>> ( istream& in, Monom& m){
    int g;
    float c;
    in >> g >> c;
    m = Monom(g, c);
    return in;
}
ostream& operator<< ( ostream& out, const Monom& m){
    if( m.coef() == 0)
        return out;
    if( m.coef() && abs(m.coef()) != 1 || m.grad() == 0) {
        if (m.coef() < 0)
            out << " - " << -m.coef();
        else
            out << " " << m.coef();
    }
    if( m.coef() == -1)
        out << "- ";
    if( m.coef() == 1)
        out << " ";
    if( abs(m.coef()) != 1 && m.grad())
        out << " * ";
    if( m.grad())
        out<< "x";
    if ( m.grad() && m.grad() != 1)
        out << "^" << m.grad() << " ";
    if ( m.grad() == 1 )
        out << " ";
    return out;
}

class Polinom{
protected:
    int nr_monoame;
    vector<Monom*> monoame;
public:
    Polinom() : nr_monoame(0){}
    Polinom(int nr_mon, const vector<Monom>& v): nr_monoame(nr_mon){
        for ( auto m: v){
            monoame.push_back(new Monom(m));
        }
        ordoneaza();
    }
    Polinom( const Polinom& P): nr_monoame(P.nr_monoame) {
        for( auto m: P.monoame){
            monoame.push_back(new Monom(*m));
        }
        ordoneaza();
    }
    virtual ~Polinom() = 0;
    int nrMonoame() const{
        return nr_monoame;
    }
    vector<Monom*> Monoame() const{
        return monoame;
    }
    void ordoneaza ();
    float valoare( float x) const;
    virtual void derivata();
    virtual float radacina( Polinom& A) const;
};

Polinom::~Polinom(){
    for(auto m: monoame){
        delete(m);
    }
    monoame.clear();
}
void Polinom:: ordoneaza() {
    Monom* v[MAX] = {};
    for( auto m: monoame){
        v[m->grad()] = m;
    }
    monoame.clear();
    for(auto & i : v)
        if( i != nullptr)
            monoame.push_back(i);
}
void Polinom:: derivata(){
    if( monoame[0]->grad() == 0){
        delete monoame[0];
        for( int i = 1; i < nr_monoame; i++ )
            monoame[i-1] = monoame[i];
        nr_monoame --;
        monoame.pop_back();
    }
    for( auto m: monoame) {
        m -> setCoef(m->coef() * m->grad());
        m -> setGrad(m->grad()-1);
    }
}
float Polinom::valoare( float x) const{
    float v = 0;
    for( auto m: monoame ){
        v += m->coef() * float(pow(x, m->grad()));
    }
    return v;
}
float Polinom::radacina( Polinom& deriv) const{
    float x = 0;
    switch (Monoame()[nrMonoame()-1]->coef() > 0){
        case false: {
            while (valoare(x + 0.01) >= 0)
                x += 0.01;
            while (valoare(x - 0.01) <= 0)
                x -= 0.01;
            break;
        }
        default: {
            while (valoare(x + 0.01) <= 0)
                x += 0.01;
            while (valoare(x - 0.01) > 0)
                x -= 0.01;
        }
    }
    int i = 0;
    while( valoare(x) != 0 && i < MAX_ITER && deriv.valoare(x) != 0) {
        x = x - valoare(x) / deriv.valoare(x);
        i++;
    }
    return x;
}

class Polinom_ireductibil : public Polinom{
public:
    Polinom_ireductibil(): Polinom(){}
    Polinom_ireductibil( int gr, const vector<Monom>& v) : Polinom(gr, v){
        try{
            test_ireductibilitate();
        }catch(int){
            throw exception();
        }
    }
    Polinom_ireductibil( const  Polinom_ireductibil& Pi ) = default;
    ~Polinom_ireductibil() override {}
    void test_ireductibilitate();
    Polinom_ireductibil& operator= ( const Polinom_ireductibil& Pi);
    bool Criteriu_Eisenstein ();
    friend istream& operator>>(istream &in, Polinom_ireductibil& A);
};

void Polinom_ireductibil::test_ireductibilitate() {
    try {
        if(!Criteriu_Eisenstein())
            throw 99;
    }
    catch( ExceptieEisenstein& e){
        cout << e.what();
        throw 99;
    }catch (int){
        if (monoame[nr_monoame-1]->grad() % 2 == 1) throw 20;
        int semn = monoame[nr_monoame-1]->coef() / abs(monoame[nr_monoame-1]->coef());
        float x = 0;
        if (semn == -1) {
            int i = 0;
            while ( valoare( x + 0.1) > valoare(x) || valoare(x) < 0) {
                x += 0.1;
                i++;
                if ( i > 1000)
                    break;
            }
            i = 0;
            while( valoare( x - 0.1) > valoare(x) || valoare(x) < 0) {
                x -= 0.1;
                i++;
                if ( i > 1000)
                    break;
            }
        } else {
            int i = 0;
            while( valoare( x - 0.1) < valoare(x) || valoare(x) > 0) {
                x -= 0.1;
                i++;
                if ( i > 1000)
                    break;
            }
            i = 0;
            while( valoare( x + 0.1) < valoare(x) || valoare(x) > 0) {
                x += 0.1;
                i++;
                if ( i > 1000)
                    break;
            }
        }
        //cout << "\nx = " << x << "; f(x) = " << valoare(x) << "\n";
        if (valoare(x) * semn <= 0)  throw 20 ;
    }
}
bool Polinom_ireductibil::Criteriu_Eisenstein() {
    int zecimale = 0, c = 0, prim = 2;
    for( auto m: monoame ){
        float aux = m->coef();
        c = 0;
        while ( aux * 10 <= INT32_MAX && aux - int(aux) != 0){
            c++;
            aux *= 10;
        }
        if( c > zecimale )
            zecimale = c;
        if(aux * 10 >= INT_MAX) {
            throw excpE;
        }
    }
    bool ok = false;
    while( prim < MAX_CIUR && !ok){
        if( ciur[prim] == 0) {
            ok = true;
            for (int m = 0; m < nr_monoame - 1; m++) {
                if (int(monoame[m]->coef() * pow(10, zecimale)) % prim != 0) {
                    ok = false;
                    break;
                }
            }
        }
        prim++;
    }
    return int(monoame[nr_monoame-1]->coef() * pow(10, zecimale)) % (prim-1) && ok && (int(monoame[0]->coef() * pow(10, zecimale)) % ((prim-1) * (prim-1))) != 0;
}
Polinom_ireductibil& Polinom_ireductibil::operator= ( const Polinom_ireductibil& Pi){
    Polinom::operator=(Pi);
    return *this;
}
ostream& operator<<( ostream& out, const Polinom_ireductibil& Pi){
    out << "Polinom ireductibil de gradul " << Pi.Monoame()[Pi.nrMonoame()-1]->grad() << ":";
    int i = 0;
    for( auto m : Pi.Monoame() ){
        if( i && m->coef() > 0)
            out << "+";
        else
            i = 1;
        out << *m;
    }
    out <<"\n";
    return out;
}
istream& operator>>(istream &in, Polinom_ireductibil& A) {
    Monom m;
    vector<Monom> v;
    in >> A.nr_monoame;
    for( int i = 0; i < A.nr_monoame; i++ ){
        in >> m;
        A.monoame.push_back(new Monom(m));
        v.push_back(m);
    }
    A.ordoneaza();
    try{
        A.test_ireductibilitate();
    } catch (int){
        cout << "Polinomul este reductibil\n";
        abort();
    }
    return in;
}


class Polinom_reductibil: public Polinom{
public:
    Polinom_reductibil(): Polinom(){}
    Polinom_reductibil(int gr, vector<Monom> v) : Polinom(gr, v){
        try{
            Polinom_ireductibil Pr(gr, v);
            cout << "Polinom ireductibil";
            abort();
        } catch(exception){
            ordoneaza();
        }
    }
    Polinom_reductibil( const Polinom_reductibil& Pr) = default ;
    ~Polinom_reductibil() override{}
    Polinom_reductibil& operator=( const Polinom_reductibil&);
    friend istream& operator>>( istream& in, Polinom_reductibil&);
};
Polinom_reductibil& Polinom_reductibil:: operator=(const Polinom_reductibil & A) {
    Polinom::operator=(A);
    return *this;
}
istream& operator>>( istream& in, Polinom_reductibil& Pr){
    Monom m;
    vector<Monom> v;
    in >> Pr.nr_monoame;
    for( int i = 0; i < Pr.nr_monoame; i++ ){
        in >> m;
        Pr.monoame.push_back(new Monom(m));
        v.push_back(m);
    }
    try{
        Pr.ordoneaza();
        Polinom_ireductibil Pi(Pr.nrMonoame(), v);
        cout << "Polinom ireductibil citit intr-un obiect de clasa Polinom_reductibil\n";
        abort();
    } catch(exception&){}
    return in;
}
ostream& operator<<( ostream& out, const Polinom_reductibil& Pr){
    out << "Polinom reductibil de gradul " << Pr.Monoame()[Pr.nrMonoame()-1]->grad() << ": ";
    int i = 0;
    for( auto m : Pr.Monoame() ){
        if( i && m->coef() > 0 )
            out << "+";
        else
            i = 1;
        out << *m;
    }
    out <<" = ";
    Polinom_reductibil d(Pr);
    d.derivata();
    float x = - Pr.radacina(d);
    if ( x == 0)
        out << "x * (";
    else {
        out << " ( x ";
        if (x > 0)
            out << "+ " << x;
        else
            out << "- " << -x;
        out << " ) * (";
    }
    float coef;
    out << Monom(Pr.Monoame()[Pr.nrMonoame() - 1]->grad() - 1, Pr.Monoame()[Pr.nrMonoame()-1]->coef());
    coef = Pr.Monoame()[Pr.nrMonoame()-1]->coef();
    int p = Pr.nrMonoame() - 2;
    for( int i = Pr.Monoame()[Pr.nrMonoame()-1]->grad() - 2; i >= 0; i--){
        if( Pr.Monoame()[p]->grad() != i + 1 ){
            coef = - x * coef;
        } else
        {
            coef = Pr.Monoame()[p]-> coef() - x * coef;
            p--;
        }
        if( coef > 0)
            out << "+";
        out << Monom( i, coef);
    }
    out << ")";
    return out;
}

bool esteReductibil ( int grad, vector<Monom> v){
    try{
        Polinom_ireductibil a(4, v);
        return false;
    } catch( exception ){
        return true;
    }
}
bool esteIreductibil ( int grad, vector<Monom> v){
    try{
        Polinom_ireductibil a(4, v);
        return true;
    } catch( exception ){
        return false;
    }
}

int main(){
    ifstream read ("Polinom.in");
    ofstream print ("Polinom.out");
    ciur_eratostene();
    int x;
    read >> x;
    print << fixed << setprecision(2);
    cout << fixed << setprecision(2);
    switch ( x ){
        case 1: {
            vector <Polinom_reductibil*> v;
            int n;
            read >> n;
            for ( int i = 0; i < n ; i++)
            {
                Polinom_reductibil A;
                read >> A;
                v.push_back(new Polinom_reductibil(A));
            }
            for ( auto P : v ){
                print << *P << "\n";
            }
            for ( auto P : v)
                delete(P);
            break;
        }
        default: {
            vector <Polinom_ireductibil*> v;
            int n;
            read >> n;
            for (int i = 0; i < n; i++) {
                Polinom_ireductibil A;
                read >> A;
                v.push_back(new Polinom_ireductibil(A));
            }
            for (auto P : v) {
                print << *P << "\n";
            }
            for ( auto P : v)
                delete(P);
        }
    }
    return 0;
}

