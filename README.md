# BLM1022_SayisalAnaliz_proje-
## 📌 Sayısal Yöntemler Uygulaması

Bu proje, **C programlama dili** kullanılarak geliştirilmiş olup çeşitli **sayısal analiz yöntemlerini** destekleyen bir terminal uygulamasıdır. Yıldız Teknik Üniversitesi **Bilgisayar Mühendisliği** bölümü, **BLM1022 Sayısal Analiz** dersi dönem projesi kapsamında hazırlanmıştır.

### 👨‍💻 Geliştirici Bilgileri

- **İsim:** Yusuf İBİN  
- **Numara:** 24011074  
- **E-posta:** yusuf.ibin@std.edu.yildiz.edu.tr  

---

## 🧠 Desteklenen Sayısal Yöntemler

Aşağıdaki yöntemler desteklenmektedir:

1. **Bisection Yöntemi**  
2. **Regula Falsi Yöntemi**  
3. **Newton-Raphson Yöntemi**  
4. **Matrisin Tersini Alma**  
5. **Cholesky Yöntemi**  
6. **Gauss-Seidel Yöntemi**  
7. **Sayısal Türev**  
8. **Simpson Yöntemi (1/3 ve 3/8)**  
9. **Trapez Kuralı**  
10. **Gregory-Newton Enterpolasyonu**

---

## 🚀 Kullanım

Program çalıştırıldığında kullanıcı aşağıdaki işlemleri yapabilir:

- Fonksiyon tabanlı yöntemler için (örneğin Bisection, Newton-Raphson):
  - Polinom, logaritmik, trigonometrik, üstel vb. fonksiyonları string olarak terminalden girer.
  - Fonksiyonlar `^` (üs), `*` (çarpma), `_` (taban) gibi özel karakterlerle ifade edilir.
- Matris tabanlı yöntemler için:
  - Matris boyutu ve elemanları kullanıcıdan alınır.

### 📋 Örnek Girdiler

#### Bisection Yöntemi

```
Fonksiyon: x^3 - 4*x + cos(x)  
Aralık: (1, 3)  
Hata: 0.001  
```

#### Regula Falsi

```
Fonksiyon: e^(-x) - x  
Aralık: (0, 1)  
Hata: 0.001  
```

#### Newton-Raphson

```
Fonksiyon: x^3 - x - 1  
Başlangıç değeri: 1.5  
Hata: 0.0001  
```

#### Sayısal Türev

```
Fonksiyon: e^(x) * sin(x)  
Nokta: 1.5  
H: 0.1  
```

#### Simpson Yöntemi

```
Fonksiyon: ln(x)  
Alt sınır: 1  
Üst sınır: 2  
Parça sayısı: 4  
Yöntem: 1/3 kuralı  
```

---

## 🛠 Derleme ve Çalıştırma

```bash
gcc main.c -o sayisal_program -lm
./sayisal_program
```

> Not: `-lm` matematik kütüphanesini bağlamak için gereklidir.
