#!/bin/bash

# Colores para logs
GREEN='\033[0;32m'
NC='\033[0m'

echo -e "${GREEN}🚀 Iniciando despliegue en EC2...${NC}"

# 1. Actualizar sistema e instalar dependencias básicas
echo -e "${GREEN}📦 Actualizando sistema instalando dependencias...${NC}"
sudo apt-get update
sudo apt-get install -y git python3-pip python3-venv nginx curl

# 2. Instalar Node.js 20 (LTS)
echo -e "${GREEN}📦 Instalando Node.js...${NC}"
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt-get install -y nodejs

# 3. Clonar repositorio
echo -e "${GREEN}📂 Clonando repositorio...${NC}"
if [ -d "BioInformaticaProyecto1-" ]; then
    echo "El directorio ya existe. Actualizando..."
    cd BioInformaticaProyecto1-
    git pull
else
    git clone https://github.com/milith0kun/BioInformaticaProyecto1-.git
    cd BioInformaticaProyecto1-
fi

# 4. Configurar Backend
echo -e "${GREEN}🐍 Configurando Backend...${NC}"
cd backend
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install gunicorn  # Servidor de producción

# Crear .env si no existe
if [ ! -f .env ]; then
    echo "Creando .env básico..."
    cp .env.example .env
    # NOTA: Aquí deberías editar el .env real con tu API Key más tarde
fi

# Crear servicio Systemd para Backend
echo -e "${GREEN}⚙️ Creando servicio de Backend...${NC}"
cat <<EOF | sudo tee /etc/systemd/system/bio_backend.service
[Unit]
Description=Gunicorn instance to serve BioInformatics API
After=network.target

[Service]
User=ubuntu
Group=www-data
WorkingDirectory=/home/ubuntu/BioInformaticaProyecto1-/backend
Environment="PATH=/home/ubuntu/BioInformaticaProyecto1-/backend/venv/bin"
ExecStart=/home/ubuntu/BioInformaticaProyecto1-/backend/venv/bin/gunicorn -w 4 -k uvicorn.workers.UvicornWorker app.main:app --bind 0.0.0.0:8000

[Install]
WantedBy=multi-user.target
EOF

sudo systemctl daemon-reload
sudo systemctl start bio_backend
sudo systemctl enable bio_backend

# 5. Configurar Frontend
echo -e "${GREEN}⚛️ Construyendo Frontend...${NC}"
cd ../frontend
npm install
npm run build

# 6. Configurar Nginx
echo -e "${GREEN}🌐 Configurando Nginx...${NC}"
cat <<EOF | sudo tee /etc/nginx/sites-available/bioinformatica
server {
    listen 80;
    server_name _;

    # Frontend (Archivos estáticos)
    location / {
        root /home/ubuntu/BioInformaticaProyecto1-/frontend/dist;
        index index.html;
        try_files \$uri \$uri/ /index.html;
    }

    # Backend API Proxy
    location /api {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
    }
}
EOF

# Activar sitio
sudo rm -f /etc/nginx/sites-enabled/default
sudo ln -sf /etc/nginx/sites-available/bioinformatica /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl restart nginx

echo -e "${GREEN}✅ ¡Despliegue completado!${NC}"
echo -e "Visita tu IP pública para ver la aplicación."
